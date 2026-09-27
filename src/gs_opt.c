/* ============================================================================
* Stomatal optimisation ("gs_opt"): Sperry et al. (2017) profit maximisation,
* as implemented in De Kauwe et al. (2022), for the sunlit and shaded big
* leaves of the sub-daily model with the SPA soil (water_balance =
* hydraulics). Replaces the Emax scheme.
*
* At each iteration of the leaf temperature loop the leaf chooses the Ci (and
* so gs, E and leaf water potential) that maximises
*
*     profit = A / max(A) - (max(k) - k) / (max(k) - kcrit)
*
* over the feasible states (xylem conductance k above kcrit, gs below the
* optional cap), where k is the whole-plant xylem conductance at the leaf
* water potential, from a cumulative Weibull vulnerability curve, and the
* leaf water potential is found from E by integrating the curve from the
* root zone water potential (Newton-Raphson on the incomplete gamma
* function). Two searches over Ci in [gamma*, Cs]: a flat grid (default,
* gs_opt_n_sample points) or a golden section search after a coarse prescan
* (falling back to the flat grid when the feasible region is too narrow for
* the prescan to resolve).
*
* This follows the JULES gs_opt_dev implementation (stom_opt_jls_mod,
* xylem_hydraulics_cumulative_weibull_jls_mod) as closely as GDAY's leaf
* model allows: A is Farquhar min(Ac, Aj) - Rd with GDAY's leaf parameters,
* the CO2 at the leaf surface (Cs) is used in place of Ca because GDAY
* solves the boundary layer separately, and by default the transpiration
* costed is the one the leaf energy balance delivers, E = gv D / P with gv
* = gs in series with the leaf boundary layer and D the leaf to air vapour
* pressure deficit (gs_opt_e = total). gs_opt_e = stomatal uses gs alone,
* as JULES, whose optimiser leaves out the boundary layer and aerodynamic
* resistances that its surface energy balance then applies.
*
* The plant conductance is per unit leaf area (kp, as JULES kmax_pft), all
* leaves draw on it in parallel, so each big leaf gets kp x its leaf area.
*
* References:
* -----------
* * Sperry JS et al. (2017) Plant, Cell & Environment 40: 816–830.
* * De Kauwe MG, Sabot MEB, Medlyn BE et al. (2022) Towards species-level
*   forecasts of drought-induced tree mortality risk. New Phytologist 235:
*   94-110, https://doi.org/10.1111/nph.18129
*
* =========================================================================== */
#include "gs_opt.h"

#define GS_OPT_GOLDEN_RATIO 0.6180339887498949
#define GS_OPT_MIN_GOOD_PRESCAN 2
#define GS_OPT_CLOSED_GSC 1E-09

typedef struct {
    /* big leaf photosynthesis (umol m-2 s-1) */
    double vcmax, Vj, km, gamma_star, rd;
    double cs;         /* CO2 at the leaf surface (umol mol-1) */
    double dleaf;      /* leaf to air VPD (Pa) */
    double gbv;        /* big leaf boundary layer conductance for H2O
                          (mol m-2 s-1), <= 0 ignores it (JULES) */
    double press;      /* Pa */
    double psi_rz;     /* root zone water potential (MPa) */
    double kmax;       /* big leaf whole plant xylem conductance (mmol m-2
                          s-1 MPa-1) */
    double kcrit;
    int    nseg;       /* 1, or root, stem & leaf in series */
    double kseg[N_PLANT_SEG];  /* maximum conductance of each segment */
    double bseg[N_PLANT_SEG];  /* Weibull parameters (MPa, -) */
    double cseg[N_PLANT_SEG];
    double gsw_max;    /* big leaf cap on gs for H2O (mol m-2 s-1), <= 0 off */
    const gs_opt_canopy_in *canopy;  /* A(Ci) & E(gs) supplied (daily), or
                                        NULL for the Farquhar leaf above */
} leaf_in;

typedef struct {
    double ci, an, gsc, e, psi, kl, psi_stem;
    int    feasible;
} leaf_state;


static void weibull_bc(double p50, double p88, double *b, double *c) {
    /*
        Cumulative Weibull k(psi) = kmax exp(-(psi / b)^c) through P50 and
        P88, as JULES pftparm_io.
    */
    *c = log(log(1.0 - 0.5) / log(1.0 - 0.88)) / log(p50 / p88);
    *b = p50 / pow(-log(1.0 - 0.5), 1.0 / *c);

    return;
}

void weibull_params(params *p, double *b, double *c) {
    weibull_bc(p->p50, p->p88, b, c);
}

static double weibull_k(double psi, double kmax, double b, double c) {
    return (kmax * exp(-pow(MAX(0.0, psi / b), c)));
}

static double lower_incomplete_gamma(double a, double x) {
    /*
        gamma(a, x) = integral(t^(a-1) e^-t dt, 0, x), series representation
        (Numerical Recipes), as JULES incomplete_gamma.
    */
    double step, sum;
    int    i;

    if (x <= 0.0) {
        return (0.0);
    }
    step = 1.0 / a;
    sum = step;
    for (i = 1; i <= 40; i++) {
        step *= x / (a + i);
        sum += step;
        if (fabs(step) < 1E-08 * fabs(sum)) {
            break;
        }
    }

    return (sum * exp(-x) * pow(x, a));
}

static double supply_psi_leaf(double e, double psi_rz, double kmax,
                              double kcrit, double b, double c, double *kl) {
    /*
        Leaf water potential (MPa) that supplies transpiration e (mmol m-2
        s-1) from the root zone, E = integral of k(psi) from psi_leaf to
        psi_rz = kmax (-b/c) [gamma(1/c, (psi_l/b)^c) - gamma(1/c,
        (psi_r/b)^c)], solved by Newton-Raphson (JULES psi_aprox_NR, 4
        iterations). Returns the conductance at psi_leaf in kl.

        A demand the curve can never supply doesn't converge, psi_leaf keeps
        falling and kl ends below kcrit, so the state is infeasible.
    */
    double a = 1.0 / c, g_root, k, kprev, e_cur, psi, klim = 0.1 * kcrit;
    int    i;

    g_root = lower_incomplete_gamma(a, pow(MAX(0.0, psi_rz / b), c));
    k = weibull_k(psi_rz, kmax, b, c);
    if (e <= 0.0 || k <= 0.0) {
        *kl = k;
        return (psi_rz);
    }
    psi = psi_rz - e / k;

    for (i = 0; i < 4; i++) {
        kprev = k;
        k = weibull_k(psi, kmax, b, c);
        e_cur = (lower_incomplete_gamma(a, pow(MAX(0.0, psi / b), c)) -
                 g_root) * kmax * (-b / c);
        psi -= (e - e_cur) / MAX(k, 1E-30);
        psi = MAX(psi, psi_rz - 5.0 * fabs(b));
        if (fabs(k - kprev) < klim) {
            break;
        }
    }
    *kl = weibull_k(psi, kmax, b, c);

    return (psi);
}

static double plant_supply(const leaf_in *in, double e, double psi_rz,
                           double *psi_stem, double *k_plant) {
    /*
        Leaf water potential (MPa) supplying e (mmol m-2 s-1) through the
        plant's segments in series (root, stem, leaf, or the whole plant as
        one), each with its own vulnerability curve. Returns the water
        potential at the top of the stem (= leaf for a single segment) and
        the whole plant conductance k = -dE/dpsi_leaf, which the profit
        is costed on:

            E = int_{psi_out}^{psi_in} k_s  =>  dpsi_out/dE =
                (k_s(psi_in) dpsi_in/dE - 1) / k_s(psi_out)

        from dpsi_soil/dE = 0, i.e. 1 / k = sum of the segments' 1 / k_s
        with no flow, and k = k(psi_leaf) for a single segment.
    */
    double psi_in = psi_rz, psi_out, k_in, k_out, dpsi_de = 0.0;
    int    s;

    *psi_stem = psi_rz;
    for (s = 0; s < in->nseg; s++) {
        k_in = weibull_k(psi_in, in->kseg[s], in->bseg[s], in->cseg[s]);
        psi_out = supply_psi_leaf(e, psi_in, in->kseg[s],
                                  in->kcrit / in->kmax * in->kseg[s],
                                  in->bseg[s],
                                  in->cseg[s], &k_out);
        dpsi_de = (k_in * dpsi_de - 1.0) / MAX(k_out, 1E-30);
        if (s == in->nseg - 2) {
            *psi_stem = psi_out;
        }
        psi_in = psi_out;
    }
    if (in->nseg == 1) {
        *psi_stem = psi_in;
    }
    *k_plant = -1.0 / dpsi_de;

    return (psi_in);
}

static void eval_ci(const leaf_in *in, double ci, leaf_state *st) {
    /* leaf state at a given Ci */
    double Ac, Aj, gsw, gv, dcs, e;

    st->ci = ci;
    if (in->canopy != NULL) {
        const gs_opt_canopy_in *cp = in->canopy;
        st->an = cp->assim(ci, cp->ctx);
        st->gsc = st->an / MAX(in->cs - ci, 0.1);
        e = MAX(0.0, cp->trans(MAX(0.0, st->gsc), cp->ctx));
        st->e = cp->e_scale * e;
        st->psi = plant_supply(in, st->e, in->psi_rz, &st->psi_stem, &st->kl);
        st->feasible = (st->kl > in->kcrit) &&
                       (st->psi <= in->psi_rz + 1E-12) &&
                       (cp->e_max < 0.0 || e <= cp->e_max);
        return;
    }

    Ac = in->vcmax * (ci - in->gamma_star) / (ci + in->km);
    Aj = in->Vj * (ci - in->gamma_star) / (ci + 2.0 * in->gamma_star);
    st->an = MIN(Ac, Aj) - in->rd;

    // floor Cs - Ci, near Cs a narrow golden bracket would otherwise give an
    // infinite gs (JULES uses 1e-2 Pa ~ 0.1 umol mol-1)
    dcs = MAX(in->cs - ci, 0.1);
    st->gsc = st->an / dcs;                       /* mol CO2 m-2 s-1 */
    gsw = GSVGSC * st->gsc;
    // the transpiration this gs gives: gs in series with the boundary layer
    // (what the energy balance delivers), or gs alone (JULES)
    gv = gsw > 0.0 && in->gbv > 0.0 ? gsw * in->gbv / (gsw + in->gbv) : gsw;
    st->e = MAX(0.0, gv * in->dleaf / in->press) * MOL_2_MMOL;
    st->psi = plant_supply(in, st->e, in->psi_rz, &st->psi_stem, &st->kl);
    st->feasible = (st->kl > in->kcrit) && (st->psi <= in->psi_rz + 1E-12) &&
                   (in->gsw_max <= 0.0 || gsw <= in->gsw_max);

    return;
}

static double profit(const leaf_state *st, double max_an, double max_kl,
                     double kcrit) {
    return (st->an / max_an - (max_kl - st->kl) / (max_kl - kcrit));
}

static void closed_state(const leaf_in *in, leaf_state *st) {
    st->ci = in->cs;
    st->an = -in->rd;
    st->gsc = GS_OPT_CLOSED_GSC;
    st->e = 0.0;
    st->psi = plant_supply(in, 0.0, in->psi_rz, &st->psi_stem, &st->kl);
    st->feasible = FALSE;
}

static void search_flat(const leaf_in *in, int n, leaf_state *best) {
    /* flat grid over [gamma*, Cs), as JULES ci_search_flat */
    leaf_state *st;
    double lo = MAX(in->gamma_star, 0.0), hi = in->cs, max_an = -1E30;
    double max_kl = -1E30, f, best_f = -1E30;
    int    i, any = FALSE;

    if ((st = (leaf_state *)malloc(n * sizeof(leaf_state))) == NULL) {
        fprintf(stderr, "gs_opt: error allocating the Ci samples\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n; i++) {
        eval_ci(in, lo + i * (hi - lo) / n, &st[i]);
        if (st[i].feasible) {
            any = TRUE;
            max_an = MAX(max_an, st[i].an);
            max_kl = MAX(max_kl, st[i].kl);
        }
    }

    closed_state(in, best);
    // no carbon to gain from opening (e.g. dawn/dusk), stay closed
    if (any && max_an > 0.0) {
        for (i = 0; i < n; i++) {
            if (st[i].feasible) {
                f = profit(&st[i], max_an, max_kl, in->kcrit);
                if (f > best_f) {
                    best_f = f;
                    *best = st[i];
                }
            }
        }
    }
    free(st);

    return;
}

static int search_golden(const leaf_in *in, int n_prescan, int n_iter,
                         leaf_state *best) {
    /*
        Golden section search for the profit maximising Ci, after a coarse
        prescan that provides max(A) and max(k) for the normalisation, as
        JULES stom_opt_golden_search. The best state seen anywhere in the
        search is kept. Returns TRUE when the feasible region is too narrow
        for the prescan to resolve and the flat search should be used.
    */
    leaf_state st, sc, sd, *pre;
    double lo = MAX(in->gamma_star, 0.0), hi = in->cs;
    double max_an = -1E30, max_kl = -1E30, best_f = -1E30, fc, fd, f;
    double a, b, cpt, dpt;
    int    i, n_good_pos = 0, any = FALSE;

    if ((pre = (leaf_state *)malloc(n_prescan * sizeof(leaf_state))) == NULL) {
        fprintf(stderr, "gs_opt: error allocating the Ci prescan\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < n_prescan; i++) {
        eval_ci(in, lo + i * (hi - lo) / n_prescan, &pre[i]);
        if (pre[i].feasible) {
            any = TRUE;
            max_an = MAX(max_an, pre[i].an);
            max_kl = MAX(max_kl, pre[i].kl);
            if (pre[i].an > 0.0) {
                n_good_pos++;
            }
        }
    }

    // under severe stress the feasible region is a sliver just above
    // gamma*, hand it to the flat search (unless even gamma* is infeasible,
    // when nothing is)
    if (pre[0].feasible && n_good_pos < GS_OPT_MIN_GOOD_PRESCAN) {
        free(pre);
        return (TRUE);
    }

    closed_state(in, best);
    if (!any || max_an <= 0.0) {
        free(pre);
        return (FALSE);
    }
    for (i = 0; i < n_prescan; i++) {
        if (pre[i].feasible) {
            f = profit(&pre[i], max_an, max_kl, in->kcrit);
            if (f > best_f) {
                best_f = f;
                *best = pre[i];
            }
        }
    }
    free(pre);

    #define GS_OPT_F(s) ((s).feasible ? profit(&(s), max_an, max_kl, \
                                               in->kcrit) : -1E30)
    a = lo;
    b = hi;
    cpt = b - GS_OPT_GOLDEN_RATIO * (b - a);
    dpt = a + GS_OPT_GOLDEN_RATIO * (b - a);
    eval_ci(in, cpt, &sc);
    eval_ci(in, dpt, &sd);
    fc = GS_OPT_F(sc);
    fd = GS_OPT_F(sd);
    if (fc > best_f) { best_f = fc; *best = sc; }
    if (fd > best_f) { best_f = fd; *best = sd; }

    for (i = 0; i < n_iter; i++) {
        // ">=": when both probes are infeasible move towards lower Ci, where
        // the feasible region is
        if (fc >= fd) {
            b = dpt;
            dpt = cpt; sd = sc; fd = fc;
            cpt = b - GS_OPT_GOLDEN_RATIO * (b - a);
            eval_ci(in, cpt, &st);
            sc = st; fc = GS_OPT_F(sc);
            f = fc;
        } else {
            a = cpt;
            cpt = dpt; sc = sd; fc = fd;
            dpt = a + GS_OPT_GOLDEN_RATIO * (b - a);
            eval_ci(in, dpt, &st);
            sd = st; fd = GS_OPT_F(sd);
            f = fd;
        }
        if (f > best_f) {
            best_f = f;
            *best = st;
        }
    }
    #undef GS_OPT_F

    return (FALSE);
}

static void setup_plant(control *c, params *p, double kmax, leaf_in *in) {
    /*
        Whole plant, or root, stem and leaf segments sharing the whole plant
        resistance (seg_frac_*), so the segments in series give kmax when
        well watered. Segment P50/P88 default to the whole plant values.
    */
    double frac[N_PLANT_SEG] = {p->seg_frac_root, p->seg_frac_stem,
                                p->seg_frac_leaf};
    double p50[N_PLANT_SEG] = {p->p50_root, p->p50_stem, p->p50_leaf};
    double p88[N_PLANT_SEG] = {p->p88_root, p->p88_stem, p->p88_leaf};
    double sum = frac[0] + frac[1] + frac[2];
    int    s;

    in->kmax = kmax;
    in->kcrit = (1.0 - p->kcrit_frac) * kmax;
    if (c->plant_segments == N_PLANT_SEG) {
        in->nseg = N_PLANT_SEG;
        for (s = 0; s < N_PLANT_SEG; s++) {
            in->kseg[s] = kmax * sum / frac[s];
            weibull_bc(p50[s] < -900.0 ? p->p50 : p50[s],
                       p88[s] < -900.0 ? p->p88 : p88[s],
                       &in->bseg[s], &in->cseg[s]);
        }
    } else {
        in->nseg = 1;
        in->kseg[0] = kmax;
        weibull_bc(p->p50, p->p88, &in->bseg[0], &in->cseg[0]);
    }

    return;
}

static void setup_leaf(control *c, canopy_wk *cw, met *m, params *p,
                       state *s, double psi_rz, leaf_in *in) {
    double jmax, J, tk;
    int    idx = cw->ileaf;

    leaf_photo_params(c, cw, p, s, &in->gamma_star, &in->km, &in->vcmax,
                      &jmax, &J, &in->Vj, &in->rd);
    in->canopy = NULL;
    in->cs = cw->Cs;
    in->dleaf = MAX(cw->dleaf, 0.0);
    in->gbv = c->gs_opt_e == GS_OPT_E_TOTAL ? cw->gbv_leaf[idx] : -1.0;
    in->press = m->press;
    in->psi_rz = psi_rz;
    setup_plant(c, p, p->kp * cw->lai_leaf[idx], in);

    // cap on gs for H2O, m s-1 -> mol m-2 s-1, scaled to the big leaf as
    // Vcmax (JULES: som_gl_max * fpar)
    if (p->gs_opt_gl_max > 0.0) {
        tk = cw->tleaf[idx] + DEG_TO_KELVIN;
        in->gsw_max = p->gs_opt_gl_max * m->press / (RGAS * tk) *
                      cw->scalex[idx];
    } else {
        in->gsw_max = -1.0;
    }

    return;
}

static void optimise(control *c, const leaf_in *in, params *p,
                     leaf_state *best) {
    if (c->gs_opt_search == GS_OPT_GOLDEN) {
        if (search_golden(in, p->gs_opt_n_prescan, p->gs_opt_n_golden,
                          best)) {
            search_flat(in, p->gs_opt_n_sample, best);
        }
    } else {
        search_flat(in, p->gs_opt_n_sample, best);
    }
}

void gs_opt_leaf(control *c, canopy_wk *cw, met *m, params *p, state *s) {
    /*
        Profit maximising An, gs and leaf water potential of the current big
        leaf, given its temperature, Cs and leaf to air VPD.
    */
    leaf_in    in;
    leaf_state best;
    int        idx = cw->ileaf;

    setup_leaf(c, cw, m, p, s, s->weighted_swp, &in);
    if (in.kmax <= 0.0 || in.Vj <= 0.0 || in.vcmax <= 0.0) {
        closed_state(&in, &best);
    } else {
        optimise(c, &in, p, &best);
    }

    // an open optimum is clipped at zero, as JULES
    cw->an_leaf[idx] = best.feasible ? MAX(0.0, best.an) : best.an;
    cw->rd_leaf[idx] = in.rd;
    cw->gsc_leaf[idx] = MAX(GS_OPT_CLOSED_GSC, best.gsc);
    cw->lwp_leaf[idx] = best.psi;
    cw->psi_stem_leaf[idx] = best.psi_stem;
    cw->kl_leaf[idx] = in.kmax > 0.0 ? best.kl / cw->lai_leaf[idx] : 0.0;

    return;
}

int gs_opt_canopy(control *c, params *p, const gs_opt_canopy_in *cp,
                  double *ci, double *gsc, double *psi_leaf, double *e) {
    /*
        Profit maximising Ci of a canopy that supplies A(Ci) and E(gs)
        (the daily MATE model). Returns FALSE when the stomata stay shut (no
        feasible open state or no carbon to gain), with ci = gamma*.
        e (mmol m-2 s-1) is the costed transpiration (e_scale x E) and
        psi_leaf the leaf water potential supplying it.
    */
    leaf_in    in;
    leaf_state best;

    memset(&in, 0, sizeof(in));
    in.canopy = cp;
    in.gamma_star = cp->gamma_star;
    in.cs = cp->ca;
    in.psi_rz = cp->psi_rz;
    in.gsw_max = -1.0;
    setup_plant(c, p, cp->kmax, &in);
    if (in.kmax <= 0.0) {
        closed_state(&in, &best);
    } else {
        optimise(c, &in, p, &best);
    }
    if (!best.feasible || best.an <= 0.0) {
        *ci = cp->gamma_star;
        *gsc = 0.0;
        *psi_leaf = cp->psi_rz;
        *e = 0.0;
        return (FALSE);
    }
    *ci = best.ci;
    *gsc = best.gsc;
    *psi_leaf = best.psi;
    *e = best.e;

    return (TRUE);
}

double gs_opt_psi_leaf(control *c, canopy_wk *cw, params *p, state *s,
                       double e_mmol, double *kl) {
    /*
        Leaf water potential (MPa) supplying the big leaf's actual
        transpiration (mmol m-2 s-1) from the energy balance; returns the
        whole plant xylem conductance per unit leaf area in kl, and sets the
        stem water potential.
    */
    leaf_in in;
    double  psi, k;
    int     idx = cw->ileaf;

    setup_plant(c, p, p->kp * cw->lai_leaf[idx], &in);
    if (in.kmax <= 0.0) {
        *kl = p->kp;
        cw->psi_stem_leaf[idx] = s->weighted_swp;
        return (s->weighted_swp);
    }
    psi = plant_supply(&in, e_mmol, s->weighted_swp,
                       &cw->psi_stem_leaf[idx], &k);
    *kl = k / cw->lai_leaf[idx];

    return (psi);
}

double gs_opt_beta(control *c, canopy_wk *cw, met *m, params *p, state *s) {
    /*
        Water stress factor of the current big leaf for decomposition,
        allocation and output: net photosynthesis at the optimum relative to
        the optimum with wet soil (psi_root_zone = 0), at the same leaf state.
    */
    leaf_in    in;
    leaf_state wet;

    setup_leaf(c, cw, m, p, s, 0.0, &in);
    if (in.kmax <= 0.0 || in.Vj <= 0.0 || in.vcmax <= 0.0) {
        return (1.0);
    }
    optimise(c, &in, p, &wet);
    if (!wet.feasible || wet.an <= 0.0) {
        return (1.0);
    }

    return (MAX(0.0, MIN(1.0, MAX(0.0, cw->an_leaf[cw->ileaf]) / wet.an)));
}
