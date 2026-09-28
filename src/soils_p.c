/* ============================================================================
* Soil P flows, following the N flows in soils.c: litter P enters the
* structural and metabolic litter pools, moves through the CENTURY SOM pools
* (active, slow, passive) with the C, and the SOM P:C of new material rises
* with the labile P (as the N:C does with inorganic N). The mineral P is a
* labile pool in Langmuir equilibrium with a sorbed pool, a strongly sorbed
* pool, an occluded pool and parent material.
*
* Ported from GDAY-CNP (Jiang et al. 2019, https://github.com/mingkaijiang/
* GDAY-CNP), with these changes:
*   - labile and sorbed P are kept on the Langmuir isotherm exactly (the
*     total available P is updated and split by solving the isotherm),
*     rather than splitting each day's net flux by the isotherm slope, which
*     drifted off the isotherm and could make the sorbed pool negative;
*   - immobilisation is capped per SOM pool, as calculate_npools applies it,
*     and cut back if the labile P can't meet it, so the P balance closes;
*   - biochemical mineralisation is capped by what's left in the slow pool;
*   - the soil structural litter P (strpfloat) uses the soil, not surface,
*     litter P;
*   - atmospheric deposition goes to the labile pool (GDAY-CNP added it to
*     the parent material).
*
* References:
*   Wang, Y. P. et al. (2007) Global Biogeochem. Cycles, 21, GB1018.
*   Yang, X. et al. (2016) Biogeosciences, 13, 2689-2709.
*   Parton, W. J. et al. (1988, 1993) CENTURY.
*
* =========================================================================== */
#include "soils_p.h"

static double pc_limit(fluxes *, double, double, double, double);
static void   litter_p_inputs(control *, fluxes *, params *, double *,
                              double *, double *);
static void   som_p_effluxes(fluxes *, params *, state *);
static double calculate_pc_slope(params *, double, double);
static double biochemical_p_mineralisation(fluxes *, params *, state *);
static double ssorb_to_mineral_p(control *, params *, state *);


void calculate_psoil_flows(control *c, fluxes *f, params *p, state *s,
                           int doy) {
    /*
        Daily soil P update; call after calculate_nsoil_flows (needs the C
        flows into the SOM pools and today's N uptake)
    */
    int    cntrl_grazing = c->grazing;
    double faecesp, psurf, psoil, slope_a, slope_s, slope_p, arg, pc_a, pc_s,
           pc_p, imm_a, imm_s, imm_p, p_out_of_active, p_out_of_slow,
           p_out_of_passive,
           tot_in, tot_out, avail, deficit, cut, scale, lab_new, sorb_old,
           tol = 1E-10;

    if (c->grazing == 2 && p->disturbance_doy == doy+1) {
        c->grazing = TRUE;
    }

    /* litter P, split into structural and metabolic */
    litter_p_inputs(c, f, p, &faecesp, &psurf, &psoil);

    /* P leaving each organic pool, at the pool's P:C */
    som_p_effluxes(f, p, s);

    /* gross mineralisation: everything that left an organic pool */
    f->pgross = (f->p_surf_struct_to_slow + f->p_surf_struct_to_active +
                 f->p_soil_struct_to_slow + f->p_soil_struct_to_active +
                 f->p_surf_metab_to_active + f->p_soil_metab_to_active +
                 f->p_active_to_slow + f->p_active_to_passive +
                 f->p_slow_to_active + f->p_slow_to_passive +
                 f->p_passive_to_active);

    /*
    ** Litter pools. The metabolic (and, with strpfloat off, structural)
    ** pools release P to, or take P from, the labile pool to stay within
    ** their P:C limits (Parton 1989 fig 2); that is f->plittrelease.
    */
    f->plittrelease = 0.0;
    s->structsurfp += f->p_surf_struct_litter -
                      (f->p_surf_struct_to_slow + f->p_surf_struct_to_active);
    s->structsoilp += f->p_soil_struct_litter -
                      (f->p_soil_struct_to_slow + f->p_soil_struct_to_active);
    if (c->strpfloat == FALSE) {
        s->structsurfp += pc_limit(f, s->structsurf, s->structsurfp,
                                   1.0/p->structcp, 1.0/p->structcp);
        s->structsoilp += pc_limit(f, s->structsoil, s->structsoilp,
                                   1.0/p->structcp, 1.0/p->structcp);
    }
    s->metabsurfp += f->p_surf_metab_litter - f->p_surf_metab_to_active;
    s->metabsurfp += pc_limit(f, s->metabsurf, s->metabsurfp,
                              1.0/150.0, 1.0/80.0);
    s->metabsoilp += f->p_soil_metab_litter - f->p_soil_metab_to_active;
    s->metabsoilp += pc_limit(f, s->metabsoil, s->metabsoilp,
                              1.0/150.0, 1.0/80.0);

    /*
    ** Round off an emptied metabolic pool to zero; what's left goes on to
    ** the active pool with the rest of the pool's efflux (so no P is lost)
    */
    if (s->metabsurfp < tol) {
        f->p_surf_metab_to_active += s->metabsurfp;
        f->pgross += s->metabsurfp;
        s->metabsurfp = 0.0;
    }
    if (s->metabsoilp < tol) {
        f->p_soil_metab_to_active += s->metabsoilp;
        f->pgross += s->metabsoilp;
        s->metabsoilp = 0.0;
    }

    /*
    ** Immobilisation: new SOM takes P at a P:C that rises linearly with the
    ** labile P between pcmin (at pmin0) and pcmax (at pmincrit), capped per
    ** pool
    */
    slope_a = calculate_pc_slope(p, p->actpcmax, p->actpcmin);
    slope_s = calculate_pc_slope(p, p->slowpcmax, p->slowpcmin);
    slope_p = calculate_pc_slope(p, p->passpcmax, p->passpcmin);
    arg = s->inorglabp - p->pmin0 / M2_AS_HA * G_AS_TONNES;
    pc_a = MIN(p->actpcmax, p->actpcmin + slope_a * arg);
    pc_s = MIN(p->slowpcmax, p->slowpcmin + slope_s * arg);
    pc_p = MIN(p->passpcmax, p->passpcmin + slope_p * arg);
    imm_a = MAX(0.0, f->c_into_active * pc_a);
    imm_s = MAX(0.0, f->c_into_slow * pc_s);
    imm_p = MAX(0.0, f->c_into_passive * pc_p);
    f->pimmob = imm_a + imm_s + imm_p;

    /*
    ** biochemical (phosphatase) mineralisation of the slow pool; it depends
    ** on the N cost of P uptake, so needs the N cycle (off in a CP run)
    */
    f->p_slow_biochemical = c->ncycle ?
                            biochemical_p_mineralisation(f, p, s) : 0.0;

    /* mineral P transfers */
    f->p_par_to_min = p->p_rate_par_weather * s->inorgparp;
    f->p_atm_dep = p->p_atm_deposition;
    f->p_ssorb_to_min = ssorb_to_mineral_p(c, p, s);
    f->p_ssorb_to_occ = (s->inorgssorbp > 0.0) ?
                        p->rate_ssorb_occ * s->inorgssorbp : 0.0;
    f->p_min_to_ssorb = (s->inorgsorbp > 0.0) ?
                        p->rate_sorb_ssorb * s->inorgsorbp : 0.0;

    /* grazer urine P goes straight to the labile pool */
    f->purine = 0.0;
    if (c->grazing) {
        f->purine = MAX(0.0, f->peaten * p->fractosoilp - faecesp);
    }

    /*
    ** If the available (labile + sorbed) P can't cover today's demand, cut
    ** the leaching and sorption first, then the immobilisation (new SOM
    ** then forms below its minimum P:C)
    */
    f->pmineralisation = f->pgross - f->pimmob + f->plittrelease;
    avail = s->inorglabp + s->inorgsorbp;
    tot_in = f->pmineralisation + f->p_slow_biochemical + f->p_par_to_min +
             f->p_atm_dep + f->purine + f->p_ssorb_to_min;
    tot_out = f->puptake + f->ploss + f->p_min_to_ssorb;
    deficit = tot_out - tot_in - avail;
    if (deficit > 0.0) {
        cut = MIN(deficit, f->ploss);
        f->ploss -= cut;
        deficit -= cut;
        cut = MIN(deficit, f->p_min_to_ssorb);
        f->p_min_to_ssorb -= cut;
        deficit -= cut;
        if (deficit > 0.0 && f->pimmob > 0.0) {
            scale = MAX(0.0, 1.0 - deficit / f->pimmob);
            imm_a *= scale;
            imm_s *= scale;
            imm_p *= scale;
            f->pimmob = imm_a + imm_s + imm_p;
            f->pmineralisation = f->pgross - f->pimmob + f->plittrelease;
        }
        tot_in = f->pmineralisation + f->p_slow_biochemical +
                 f->p_par_to_min + f->p_atm_dep + f->purine +
                 f->p_ssorb_to_min;
        tot_out = f->puptake + f->ploss + f->p_min_to_ssorb;
    }

    /*
    ** SOM pools: the P flowing between pools was mineralised (pgross) and
    ** the new SOM immobilises P at its P:C
    */
    p_out_of_active = f->p_active_to_slow + f->p_active_to_passive;
    p_out_of_slow = f->p_slow_to_active + f->p_slow_to_passive +
                    f->p_slow_biochemical;
    p_out_of_passive = f->p_passive_to_active;
    s->activesoilp += imm_a - p_out_of_active;
    s->slowsoilp += imm_s - p_out_of_slow;
    s->passivesoilp += imm_p - p_out_of_passive;

    /*
    ** Mineral pools. Labile and sorbed P stay in Langmuir equilibrium:
    ** update their total, then split it on the isotherm. The strongly sorbed
    ** pool is fed from the sorbed pool.
    */
    sorb_old = s->inorgsorbp;
    avail = MAX(0.0, avail + tot_in - tot_out);
    lab_new = langmuir_labile_p(avail, p->smax, p->ks);
    f->p_lab_in = tot_in;
    f->p_lab_out = f->puptake + f->ploss;
    f->p_sorb_out = f->p_min_to_ssorb;
    s->inorglabp = lab_new;
    s->inorgsorbp = avail - lab_new;
    /* net labile -> sorbed exchange (negative when sorbed P desorbs) */
    f->p_sorb_in = s->inorgsorbp - sorb_old + f->p_min_to_ssorb;
    s->inorgavlp = s->inorglabp + s->inorgsorbp;

    s->inorgssorbp += f->p_min_to_ssorb - f->p_ssorb_to_occ -
                      f->p_ssorb_to_min;
    s->inorgoccp += f->p_ssorb_to_occ;
    s->inorgparp -= f->p_par_to_min;

    c->grazing = cntrl_grazing;

    return;
}

double langmuir_labile_p(double avail, double smax, double ks) {
    /*
        Labile P, L, in Langmuir equilibrium with the sorbed P, S, given
        their total A = L + S, where S = smax L / (ks + L). The positive root
        of L^2 + (ks + smax - A) L - A ks = 0, in a form that doesn't lose
        precision when b is large.
    */
    double b, disc;

    if (avail <= 0.0)
        return (0.0);
    if (smax <= 0.0)
        return (avail);

    b = ks + smax - avail;
    disc = sqrt(b * b + 4.0 * avail * ks);
    if (b >= 0.0)
        return (MIN(avail, 2.0 * avail * ks / (b + disc)));
    else
        return (MIN(avail, 0.5 * (disc - b)));
}

static void litter_p_inputs(control *c, fluxes *f, params *p,
                            double *faecesp, double *psurf, double *psoil) {
    /*
        Litter and faeces P, partitioned into the structural and metabolic
        pools. Structural P:C is either fixed (structcp) or a fraction of the
        metabolic P:C (strpfloat).
    */
    double denom;

    /* faeces P at the faeces C:P, no more than the grazed P returned */
    *faecesp = 0.0;
    if (c->grazing) {
        *faecesp = MIN(f->faecesc / p->faecescp, f->peaten * p->fractosoilp);
    }

    *psurf = f->deadleafp + f->deadbranchp + f->deadstemp + *faecesp;
    *psoil = f->deadrootp + f->deadcrootp;

    if (c->strpfloat) {
        denom = f->surf_struct_litter * p->structratp + f->surf_metab_litter;
        f->p_surf_struct_litter = float_eq(denom, 0.0) ? 0.0 :
            *psurf * f->surf_struct_litter * p->structratp / denom;

        denom = f->soil_struct_litter * p->structratp + f->soil_metab_litter;
        f->p_soil_struct_litter = float_eq(denom, 0.0) ? 0.0 :
            *psoil * f->soil_struct_litter * p->structratp / denom;
    } else {
        /* all of it goes to structural if there isn't enough */
        f->p_surf_struct_litter = MIN(*psurf,
                                      f->surf_struct_litter / p->structcp);
        f->p_soil_struct_litter = MIN(*psoil,
                                      f->soil_struct_litter / p->structcp);
    }
    f->p_surf_metab_litter = *psurf - f->p_surf_struct_litter;
    f->p_soil_metab_litter = *psoil - f->p_soil_struct_litter;

    return;
}

static void som_p_effluxes(fluxes *f, params *p, state *s) {
    /*
        P leaving the litter and SOM pools, split between destinations as
        the C is (see soils.c), at the source pool's P:C
    */
    double frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    double out, sigwt;

    /* structural */
    out = s->structsurfp * p->decayrate[0];
    sigwt = out / (p->ligshoot * 0.7 + (1.0 - p->ligshoot) * 0.55);
    f->p_surf_struct_to_slow = sigwt * p->ligshoot * 0.7;
    f->p_surf_struct_to_active = sigwt * (1.0 - p->ligshoot) * 0.55;

    out = s->structsoilp * p->decayrate[2];
    sigwt = out / (p->ligroot * 0.7 + (1.0 - p->ligroot) * 0.45);
    f->p_soil_struct_to_slow = sigwt * p->ligroot * 0.7;
    f->p_soil_struct_to_active = sigwt * (1.0 - p->ligroot) * 0.45;

    /* metabolic */
    f->p_surf_metab_to_active = s->metabsurfp * p->decayrate[1];
    f->p_soil_metab_to_active = s->metabsoilp * p->decayrate[3];

    /* active */
    out = s->activesoilp * p->decayrate[4];
    sigwt = out / (1.0 - frac_microb_resp);
    f->p_active_to_slow = sigwt * (1.0 - frac_microb_resp - 0.004);
    f->p_active_to_passive = sigwt * 0.004;

    /* slow */
    out = s->slowsoilp * p->decayrate[5];
    sigwt = out / 0.45;
    f->p_slow_to_active = sigwt * 0.42;
    f->p_slow_to_passive = sigwt * 0.03;

    /* passive */
    f->p_passive_to_active = s->passivesoilp * p->decayrate[6];

    return;
}

static double pc_limit(fluxes *f, double cpool, double ppool, double pcmin,
                       double pcmax) {
    /*
        P to add to a litter pool (negative: to release to the labile pool)
        to keep its P:C within pcmin to pcmax; the release is summed in
        f->plittrelease
    */
    double pmax = cpool * pcmax;
    double pmin = cpool * pcmin;

    if (ppool > pmax) {
        f->plittrelease += ppool - pmax;
        return (pmax - ppool);
    } else if (ppool < pmin) {
        f->plittrelease -= pmin - ppool;
        return (pmin - ppool);
    }
    return (0.0);
}

static double calculate_pc_slope(params *p, double pcmax, double pcmin) {
    /* slope of new SOM P:C vs labile P (per t P/ha) */
    return ((pcmax - pcmin) / (p->pmincrit - p->pmin0) *
            M2_AS_HA / G_AS_TONNES);
}

static double biochemical_p_mineralisation(fluxes *f, params *p, state *s) {
    /*
        Phosphatase mineralisation of the slow SOM pool (Wang et al. 2007),
        which starts once the N cost of P uptake (g N / g P) exceeds a
        critical value and growth is more limited by P than N. Capped at the
        P left in the slow pool after its decay.
    */
    double c_gain_of_p, c_gain_of_n, n_cost_of_p, x, rate, left;

    c_gain_of_p = (f->puptake > 0.0) ? f->npp / f->puptake : 0.0;
    c_gain_of_n = (f->nuptake > 0.0) ? f->npp / f->nuptake : 0.0;
    n_cost_of_p = (c_gain_of_n > 0.0) ? c_gain_of_p / c_gain_of_n : 0.0;

    rate = 0.0;
    if (c_gain_of_p > c_gain_of_n && n_cost_of_p > p->crit_n_cost_of_p) {
        x = n_cost_of_p - p->crit_n_cost_of_p;
        rate = p->max_p_biochemical * x / (x + p->biochemical_p_constant);
    }
    left = s->slowsoilp - (f->p_slow_to_active + f->p_slow_to_passive);

    return (MAX(0.0, MIN(rate, left)));
}

static double ssorb_to_mineral_p(control *c, params *p, state *s) {
    /*
        Strongly sorbed -> labile P: either a constant rate (psecmnp) or the
        CENTURY rate, which rises with soil pH between phmin and phmax and
        with sand content
    */
    double slope, rate;

    if (s->inorgssorbp <= 0.0)
        return (0.0);

    if (c->text_effect_p) {
        slope = (p->phtextmax - p->phtextmin) / (p->phmax - p->phmin);
        if (p->soilph < p->phmin)
            rate = p->phtextmin;
        else if (p->soilph > p->phmax)
            rate = p->phtextmax;
        else
            rate = p->phtextmin + slope * (p->soilph - p->phmin);
        rate += p->phtextslope * (1.0 - p->finesoil);
    } else {
        rate = p->psecmnp;
    }
    return (MAX(0.0, rate * s->inorgssorbp));
}
