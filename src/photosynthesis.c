/* ============================================================================
* Photosynthesis - C3 & C4
*
* see below
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   01.09.2018
*
* =========================================================================== */
#include "photosynthesis.h"
#include "gs_opt.h"
#include "water_balance.h"
#include "water_balance_sub_daily.h"

/* one MATE half day (am or pm) for the daily gs_opt */
typedef struct {
    params *p;
    state  *s;
    double  gamma_star, km, vcmax, jmax;
    double  par;       /* incident PAR (umol m-2 d-1) */
    double  apar;      /* absorbed PAR (umol m-2 d-1) */
    double  daylen;    /* h */
    double  secs;      /* seconds in the half day */
    double  press, vpd, tair, wind, rnet, ca;
} mate_half;

static double mate_half_gpp(double ci, void *ctx) {
    /* canopy GPP at Ci, mean over the half day (umol m-2 s-1) */
    mate_half *h = (mate_half *)ctx;
    double alpha, ac, aj;

    alpha = calculate_quantum_efficiency(h->p, ci, h->gamma_star);
    ac = assim(ci, h->gamma_star, h->vcmax, h->km);
    aj = assim(ci, h->gamma_star, h->jmax / 4.0, 2.0 * h->gamma_star);

    return (h->apar / 2.0 *
            epsilon(h->p, MIN(ac, aj), h->par, alpha, h->daylen) / h->secs);
}

static double mate_half_trans(double gsc, void *ctx) {
    /* canopy transpiration (mmol m-2 s-1) at gs, as the water balance */
    mate_half *h = (mate_half *)ctx;
    double ga, gsv, trans, LE, omega;

    penman_canopy_wrapper(h->p, h->s, h->press, h->vpd, h->tair, h->wind,
                          h->rnet, h->ca, 0.0, gsc, &ga, &gsv, &trans, &LE,
                          &omega);
    return (trans * MOL_2_MMOL);
}

static double mate_gs_opt(control *c, params *p, state *s, mate_half *h,
                          double psi_rz, double *gsc, double *psi_leaf) {
    /*
        Profit maximising Ci for the half day (Sperry et al. 2017; De Kauwe
        et al. 2022), with the same hydraulics as the sub-daily model: the
        canopy's plant conductance is kp x LAI, and the hydraulic cost is on
        the peak (midday) transpiration, gs_opt_peak_e x the half-day mean,
        since the xylem sees the peak and the vulnerability curve is
        nonlinear. The mean transpiration can't exceed what the bucket
        holds: half of the root zone's extractable water in each half day,
        so the GPP chosen is one the water balance can deliver (the bucket
        has no soil to root resistance to limit supply as the soil dries).
    */
    gs_opt_canopy_in cp;
    double ci, e;

    if (h->apar <= 0.0 || s->lai <= 0.0) {
        *gsc = 0.0;
        *psi_leaf = psi_rz;
        return (h->gamma_star);
    }
    cp.assim = mate_half_gpp;
    cp.trans = mate_half_trans;
    cp.ctx = h;
    cp.gamma_star = h->gamma_star;
    cp.ca = h->ca;
    cp.psi_rz = psi_rz;
    cp.kmax = p->kp * s->lai;
    cp.e_scale = p->gs_opt_peak_e;
    // mm per half day -> mmol m-2 s-1
    cp.e_max = 0.5 * s->pawater_root /
               (MOLE_WATER_2_G_WATER * G_TO_KG * MMOL_2_MOL * h->secs);
    gs_opt_canopy(c, p, &cp, &ci, gsc, psi_leaf, &e);

    return (ci);
}

static double bucket_psi(params *p, state *s) {
    /* root zone water potential (MPa) of the bucket, from the retention
       curve (soil_hydraulics) at the root zone water content */
    double theta = p->theta_wp_root + s->pawater_root / p->rooting_depth;

    return (MIN(0.0, soil_psi_raw(p, 0, theta)));
}

void leaf_photo_params(control *c, canopy_wk *cw, params *p, state *s,
                       double *gamma_star, double *km, double *vcmax,
                       double *jmax, double *J, double *Vj, double *rd) {
    //
    //  Photosynthetic parameters of the current big leaf (sunlit or shaded)
    //  at its leaf temperature, scaled from the single leaf to the big leaf
    //  (Wang & Leuning 1998 appendix C), and the electron transport rate for
    //  the leaf's absorbed PAR. Shared by the Medlyn and gs_opt models.
    //
    double tleaf = cw->tleaf[cw->ileaf];
    double scalex = cw->scalex[cw->ileaf];

    *gamma_star = calc_co2_compensation_point(p, tleaf);
    *km = calculate_michaelis_menten(p, tleaf);
    calculate_jmaxt_vcmaxt(c, cw, p, s, tleaf, jmax, vcmax);

    // leaf respiration in the light, Collatz et al. 1991
    *rd = 0.015 * *vcmax;

    // Scaling from single leaf to canopy, see Wang & Leuning 1998 appendix C
    *vcmax *= scalex;
    *jmax *= scalex;
    *rd *= scalex;

    // Rate of electron transport, which is a function of absorbed PAR
    calc_electron_transport_rate(p, cw->apar_leaf[cw->ileaf], *jmax, J, Vj);

    return;
}

void photosynthesis_C3(control *c, canopy_wk *cw, met *m, params *p, state *s) {
    //
    //  Calculate photosynthesis following Farquhar & von Caemmerer, this is an
    //  implementation of the routinue in MAESTRA
    //
    //  References:
    //  -----------
    //  * GD Farquhar, S Von Caemmerer (1982) Modelling of photosynthetic
    //    response to environmental conditions. Encyclopedia of plant
    //    physiology 12, 549-587.
    //  * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    //  * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.
    //

    double gamma_star, km, jmax, vcmax, rd, J, Vj, gs_over_a, g0;
    double Ci, Ac, Aj, Cs, dleaf, dleaf_kpa;
    int    idx, error = FALSE;
    double g0_zero = 1E-09; // numerical issues, don't use zero

    // unpack some stuff
    idx = cw->ileaf;
    Cs = cw->Cs;
    dleaf = cw->dleaf;

    leaf_photo_params(c, cw, p, s, &gamma_star, &km, &vcmax, &jmax, &J, &Vj,
                      &rd);

    // Deal with extreme cases
    if (jmax <= 0.0 || vcmax <= 0.0 || isnan(J)) {
        cw->an_leaf[idx] = -rd;
        cw->rd_leaf[idx] = rd;
        cw->gsc_leaf[idx] = g0_zero;
    } else {
        // Hardwiring this for Medlyn gs model for the moment, till I figure
        // out the best structure

        // For the medlyn model this is already in conductance to CO2, so the
        // 1.6 from the corrigendum to Medlyn et al 2011 is missing here
        dleaf_kpa = dleaf * PA_2_KPA;
        if (dleaf_kpa < 0.05) {
            dleaf_kpa = 0.05;
        }

        // This is calculated by SPA hydraulics so we don't need to account for
        // water stress on gs.
        if (c->water_balance == HYDRAULICS) {
            gs_over_a = (1.0 + p->g1 / sqrt(dleaf_kpa)) / Cs;
        } else {
            gs_over_a = (1.0 + (p->g1 * s->wtfac_root) / sqrt(dleaf_kpa)) / Cs;
        }
        g0 = g0_zero;

        // Solution when Rubisco activity is limiting
        error = solve_ci(g0, gs_over_a, rd, Cs, gamma_star, vcmax, km, &Ci);

        if (error || Ci <= 0.0 || Ci > Cs) {
            Ac = 0.0;
        } else {
            Ac = vcmax * (Ci - gamma_star) / (Ci + km);
        }

        // Solution when electron transport rate is limiting
        error = solve_ci(g0, gs_over_a, rd, Cs, gamma_star, Vj,
                         2.0*gamma_star, &Ci);

        if (error) {
            // no real root, handled by the compensation point case below
            Aj = 0.0;
        } else {
            Aj = Vj * (Ci - gamma_star) / (Ci + 2.0 * gamma_star);
        }

        // Below light compensation point?
        if (Aj - rd < 1E-6) {
            Ci = Cs;
            Aj = Vj * (Ci - gamma_star) / (Ci + 2.0 * gamma_star);
        }
        cw->an_leaf[idx] = MIN(Ac, Aj) - rd;
        cw->rd_leaf[idx] = rd;
        cw->gsc_leaf[idx] = MAX(g0, g0 + gs_over_a * cw->an_leaf[idx]);
    }

    return;
}

int calc_electron_transport_rate(params *p, double par, double jmax, double *J,
                                 double *Vj) {
    //
    //  Electron transport rate for a given absorbed irradiance
    //
    //  Reference:
    //  ----------
    //  * Farquhar G.D. & Wong S.C. (1984) An empirical model of stomatal
    //    conductance. Australian Journal of Plant Physiology 11, 191-210,
    //    eqn A but probably clearer in:
    //  * Leuning, R. et a., Leaf nitrogen, photosynthesis, conductance and
    //    transpiration: scaling from leaves to canopies, Plant Cell Environ.,
    //    18, 1183– 1200, 1995. Leuning 1995, eqn C3.
    //
    int large_root = FALSE, error = FALSE;
    double A, B, C;

    A = p->theta;
    B = -(p->alpha_j * par + jmax);
    C = p->alpha_j * par * jmax;

    *J = quad(A, B, C, large_root, &error);
    *Vj = *J / 4.0;

    return error;
}

int solve_ci(double g0, double gs_over_a, double rd, double Cs,
             double gamma_star, double gamma, double beta, double *Ci) {
    //
    //  Solve intercellular CO2 concentration using quadric equation, following
    //  Leuning 1990, see eqn 15a-c, solving simultaneous solution for Eqs 2,
    //  12 and 13
    //
    //  Reference:
    //  ----------
    //  * Leuning (1990) Modelling Stomatal Behaviour and Photosynthesis of
    //    Eucalyptus grandis. Aust. J. Plant Physiol., 17, 159-75.
    //
    int large_root = TRUE, error = FALSE;
    double A, B, C, arg1, arg2, arg3;

    A = g0 + gs_over_a * (gamma - rd);

    arg1 = (1. - Cs * gs_over_a) * (gamma - rd);
    arg2 = g0 * (beta - Cs);
    arg3 = gs_over_a * (gamma * gamma_star + beta * rd);
    B = arg1 + arg2 - arg3;

    arg1 = -(1.0 - Cs * gs_over_a);
    arg2 = (gamma * gamma_star + beta * rd);
    arg3 =  g0 * beta * Cs;
    C = arg1 * arg2 - arg3;

    *Ci = quad(A, B, C, large_root, &error);

    return (error);
}

double calc_co2_compensation_point(params *p, double tleaf) {
    //
    //  CO2 compensation point in the absence of non-photorespiratory
    //  respiration.
    //
    //  Parameters:
    //  ----------
    //  tleaf : float
    //      leaf temperature (deg C)
    //
    //  Returns:
    //  -------
    //  gamma_star : float
    //      CO2 compensation point in the abscence of mitochondrial respiration
    //      (umol mol-1)
    //
    //  References:
    //  -----------
    //  * Bernacchi et al 2001 PCE 24: 253-260
    //
    return (arrhenius(p->gamstar25, p->eag, tleaf, p->measurement_temp));
}

double calculate_michaelis_menten(params *p, double tleaf) {
    //
    //  Effective Michaelis-Menten coefficent of Rubisco activity
    //
    //  Parameters:
    //  ----------
    //  tleaf : float
    //      leaf temperature (deg C)
    //
    //  Returns:
    //  -------
    //  Km : float
    //      Effective Michaelis-Menten constant for Rubisco catalytic activity
    //      (umol mol-1)
    //
    //  References:
    //  -----------
    //  Rubisco kinetic parameter values are from:
    //  * Bernacchi et al. (2001) PCE, 24, 253-259.
    //  * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    double Kc, Ko, Km;

    // Michaelis-Menten coefficents for carboxylation by Rubisco
    Kc = arrhenius(p->kc25, p->eac, tleaf, p->measurement_temp);

    // Michaelis-Menten coefficents for oxygenation by Rubisco
    Ko = arrhenius(p->ko25, p->eao, tleaf, p->measurement_temp);

    // return effective Michaelis-Menten coefficient for CO2
    Km = Kc * (1.0 + p->oi / Ko);

    return (Km);
}

void calculate_jmaxt_vcmaxt(control *c, canopy_wk *cw, params *p, state *s,
                            double tleaf, double *jmax, double *vcmax) {
    //
    //  Calculate the potential electron transport rate (Jmax) and the
    //  maximum Rubisco activity (Vcmax) at the leaf temperature.
    //
    //  For Jmax -> peaked arrhenius is well behaved for tleaf < 0.0
    //
    //  Parameters:
    //  ----------
    //  tleaf : float
    //      air temperature (deg C)
    //  jmax : float
    //      the potential electron transport rate at the leaf temperature
    //      (umol m-2 s-1)
    //  vcmax : float
    //      the maximum Rubisco activity at the leaf temperature (umol m-2 s-1)
    //
    double jmax25, vcmax25;
    double lower_bound = p->photo_tlow;
    double upper_bound = p->photo_thigh;
    double tref = p->measurement_temp;

    // NB. the same top-of-canopy values are used for both leaves, the
    // sunlit/shaded difference comes from scalex in photosynthesis_C3
    if (c->modeljm == 0) {
        *jmax = p->jmax;
        *vcmax = p->vcmax;
    } else if (c->modeljm == 1) {
        vcmax25 = (p->vcmaxna * cw->N0 + p->vcmaxnb);
        jmax25 = (p->jmaxna * cw->N0 + p->jmaxnb);
        *vcmax = vcmax_temperature(p, vcmax25, tleaf);
        *jmax = peaked_arrhenius(jmax25, p->eaj, tleaf, tref, p->delsj, p->edj);
    } else if (c->modeljm == 2) {
        vcmax25 = (p->vcmaxna * cw->N0 + p->vcmaxnb);
        jmax25 = (p->jv_slope * vcmax25 - p->jv_intercept);
        *vcmax = vcmax_temperature(p, vcmax25, tleaf);
        *jmax = peaked_arrhenius(jmax25, p->eaj, tleaf, tref, p->delsj, p->edj);
    } else if (c->modeljm == 3) {
        jmax25 = p->jmax;
        vcmax25 = p->vcmax;
        *vcmax = vcmax_temperature(p, vcmax25, tleaf);
        *jmax = peaked_arrhenius(jmax25, p->eaj, tleaf, tref, p->delsj, p->edj);
    } else {
        fprintf(stderr, "You haven't set Jmax/Vcmax model: modeljm \n");
        exit(EXIT_FAILURE);
    }

    // Bucket model water stress acts through g1 (photosynthesis_C3) and,
    // optionally, on Jmax/Vcmax as a proxy for non-stomatal limitation
    if (c->water_balance == BUCKET && c->nonstomatal_limitation) {
        *jmax *= s->wtfac_root;
        *vcmax *= s->wtfac_root;
    }

    // Jmax/Vcmax forced linearly to zero at low T
    if (tleaf < lower_bound) {
        *jmax = 0.0;
        *vcmax = 0.0;
    } else if (tleaf < upper_bound) {
        *jmax *= (tleaf - lower_bound) / (upper_bound - lower_bound);
        *vcmax *= (tleaf - lower_bound) / (upper_bound - lower_bound);
    }

    return;
}

double calc_leaf_day_respiration(double tleaf, double Rd0) {
    //
    //  Calculate leaf respiration in the light using a Q10 (exponential)
    //  formulation
    //
    //  Parameters:
    //  ----------
    //  Rd0 : float
    //      Estimate of respiration rate at the reference temperature 25 deg C
    //      or or 298 K
    //  Tleaf : float
    //      leaf temperature
    //
    //  Returns:
    //  -------
    //  Rd : float
    //      leaf respiration in the light
    //

    double Rd;
    // tbelow is the minimum temperature at which respiration occurs;
    //   below this, Rd(T) is set to zero.
    double tbelow = -100.0;

    // Determines the respiration in the light, relative to
    //   respiration in the dark. For example, if DAYRESP = 0.7, respiration
    //   in the light is 30% lower than in the dark.
    double day_resp = 1.0;  // no effect of light

    // Default temp at which maint resp specified
    double rtemp = 0.0;

    // exponential coefficient of the temperature response of foliage
    //   respiration
    double q10f = 0.067;

    if (tleaf >= tbelow)
        Rd = Rd0 * day_resp * exp(q10f * (tleaf - rtemp));
    else
        Rd = 0.0;

    return (Rd);
}

double arrhenius(double k25, double Ea, double T, double Tref) {
    //
    //  Temperature dependence of kinetic parameters is described by an
    //  Arrhenius function
    //
    //  Parameters:
    //  ----------
    //  k25 : float
    //      rate parameter value at 25 degC
    //  Ea : float
    //      activation energy for the parameter [J mol-1]
    //  T : float
    //      temperature [deg C]
    //  Tref : float
    //      measurement temperature [deg C]
    //
    //  Returns:
    //  -------
    //  kt : float
    //      temperature dependence on parameter
    //
    //  References:
    //  -----------
    //  * Medlyn et al. 2002, PCE, 25, 1167-1179.
    //
    double Tk, TrefK;
    Tk = T + DEG_TO_KELVIN;
    TrefK = Tref + DEG_TO_KELVIN;

    return (k25 * exp(Ea * (T - Tref) / (RGAS * Tk * TrefK)));
}

double peaked_arrhenius(double k25, double Ea, double T, double Tref,
                        double deltaS, double Hd) {
    //
    //  Temperature dependancy approximated by peaked Arrhenius eqn,
    //  accounting for the rate of inhibition at higher temperatures.
    //
    //  Parameters:
    //  ----------
    //  k25 : float
    //      rate parameter value at 25 degC
    //  Ea : float
    //      activation energy for the parameter [J mol-1]
    //  T : float
    //      temperature [deg C]
    //  Tref : float
    //      measurement temperature [deg C]
    //  deltaS : float
    //      entropy factor [J mol-1 K-1)
    //  Hd : float
    //      describes rate of decrease about the optimum temp [J mol-1]
    //
    //  Returns:
    //  -------
    //  kt : float
    //      temperature dependence on parameter
    //
    //  References:
    //  -----------
    //  * Medlyn et al. 2002, PCE, 25, 1167-1179.
    //

    double arg1, arg2, arg3;
    double Tk, TrefK;
    Tk = T + DEG_TO_KELVIN;
    TrefK = Tref + DEG_TO_KELVIN;

    arg1 = arrhenius(k25, Ea, T, Tref);
    arg2 = 1.0 + exp((deltaS * TrefK - Hd) / (RGAS * TrefK));
    arg3 = 1.0 + exp((deltaS * Tk - Hd) / (RGAS * Tk));

    return (arg1 * arg2 / arg3);
}


double vcmax_temperature(params *p, double vcmax25, double tleaf) {
    //
    //  Vcmax at leaf (or air) temperature (deg C): plain Arrhenius, or
    //  peaked Arrhenius (Medlyn et al. 2002) when a deactivation energy (edv)
    //  is given
    //
    if (p->edv > 0.0) {
        return (peaked_arrhenius(vcmax25, p->eav, tleaf, p->measurement_temp,
                                 p->delsv, p->edv));
    } else {
        return (arrhenius(vcmax25, p->eav, tleaf, p->measurement_temp));
    }
}

double quad(double a, double b, double c, bool large, int *error) {
    //
    //  quadratic solution
    //
    //  Parameters:
    //  ----------
    //  a : float
    //      co-efficient
    //  b : float
    //      co-efficient
    //  c : float
    //      co-efficient
    //
    //  Returns:
    //  -------
    //  root : float
    //

    double d, root;

    // discriminant
    d = (b * b) - 4.0 * a * c;
    if (d < 0.0) {
        //fprintf(stderr, "imaginary root found\n");
        // exit(EXIT_FAILURE);
        *error = TRUE;
    }

    if (large) {
        if (a == 0.0 && b > 0.0) {
            root = -c / b;
        } else if (a == 0.0 && b == 0.0) {
            root = 0.0;
            if (c != 0.0) {
                // fprintf(stderr, "Can't solve quadratic\n");
                // exit(EXIT_FAILURE);
                *error = TRUE;
            }
        } else {
            root = (-b + sqrt(d)) / (2.0 * a);
        }

    } else {
        if (a == 0.0 && b > 0.0) {
            root = -c / b;
        } else if (a == 0.0 && b == 0.0) {
            root = 0.0;
            if (c != 0.0) {
                // fprintf(stderr, "Can't solve quadratic\n");
                //exit(EXIT_FAILURE);
                *error = TRUE;
            }
        } else {
            root = (-b - sqrt(d)) / (2.0 * a);
        }
    }
    return (root);
}


void mate_C3_photosynthesis(control *c, fluxes *f, met *m, params *p, state *s,
                            double daylen, double ncontent) {
    /*

    MATE simulates big-leaf C3 photosynthesis (GPP) based on Sands (1995),
    accounting for diurnal variations in irradiance and temp (am [sunrise-noon],
    pm[noon to sunset]).

    References:
    -----------
    * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    * McMurtrie, R. E. et al. (2008) Functional Change Biology, 35, 521-34.
    * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.

    Rubisco kinetic parameter values are from:
    * Bernacchi et al. (2001) PCE, 24, 253-259.
    * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    */
    double N0, gamma_star_am,
           gamma_star_pm, Km_am, Km_pm, jmax_am, jmax_pm, vcmax_am, vcmax_pm,
           ci_am, ci_pm, alpha_am, alpha_pm, ac_am, ac_pm, aj_am, aj_pm,
           asat_am, asat_pm, lue_am, lue_pm, lue_avg, conv;
    double psi_rz, psi_am, psi_pm, frac_canopy;
    mate_half h;
    double mt = p->measurement_temp + DEG_TO_KELVIN;

    /* Calculate mate params & account for temperature dependencies */
    N0 = calculate_top_of_canopy_n(p, s, ncontent);

    gamma_star_am = calculate_co2_compensation_point(p, m->Tk_am, mt);
    gamma_star_pm = calculate_co2_compensation_point(p, m->Tk_pm, mt);

    Km_am = calculate_michaelis_menten_parameter(p, m->Tk_am, mt);
    Km_pm = calculate_michaelis_menten_parameter(p, m->Tk_pm, mt);

    calculate_jmax_and_vcmax(c, p, s, m->Tk_am, N0, &jmax_am, &vcmax_am, mt);
    calculate_jmax_and_vcmax(c, p, s, m->Tk_pm, N0, &jmax_pm, &vcmax_pm, mt);

    /* Covert PAR units (umol PAR MJ-1) */
    conv = MJ_TO_J * J_2_UMOL;
    m->par *= conv;

    if (c->gs_model == GS_OPT) {
        /* profit maximising Ci for each half day */
        psi_rz = bucket_psi(p, s);
        frac_canopy = 1.0 - exp(-0.398 * s->lai);   /* as the water balance */
        h.p = p;
        h.s = s;
        h.par = m->par;
        h.apar = s->lai > 0.0 ? m->par * s->fipar : 0.0;
        h.daylen = daylen;
        h.secs = SECS_IN_HOUR * daylen / 2.0;
        h.press = m->press;
        h.ca = m->Ca;

        h.gamma_star = gamma_star_am;
        h.km = Km_am;
        h.vcmax = vcmax_am;
        h.jmax = jmax_am;
        h.vpd = m->vpd_am;
        h.tair = m->tair_am;
        h.wind = m->wind_am;
        h.rnet = calc_net_radiation(c, p, m->sw_rad_am, m->tair_am,
                                    m->lwdown_am) * frac_canopy;
        ci_am = mate_gs_opt(c, p, s, &h, psi_rz, &f->gsc_am, &psi_am);

        h.gamma_star = gamma_star_pm;
        h.km = Km_pm;
        h.vcmax = vcmax_pm;
        h.jmax = jmax_pm;
        h.vpd = m->vpd_pm;
        h.tair = m->tair_pm;
        h.wind = m->wind_pm;
        h.rnet = calc_net_radiation(c, p, m->sw_rad_pm, m->tair_pm,
                                    m->lwdown_pm) * frac_canopy;
        ci_pm = mate_gs_opt(c, p, s, &h, psi_rz, &f->gsc_pm, &psi_pm);

        s->predawn_swp = psi_rz;
        s->midday_lwp = MIN(psi_am, psi_pm);
    } else {
        ci_am = calculate_ci(c, p, s, m->vpd_am, m->Ca);
        ci_pm = calculate_ci(c, p, s, m->vpd_pm, m->Ca);
    }

    /* quantum efficiency calculated for C3 plants */
    alpha_am = calculate_quantum_efficiency(p, ci_am, gamma_star_am);
    alpha_pm = calculate_quantum_efficiency(p, ci_pm, gamma_star_pm);

    /* Rubisco carboxylation limited rate of photosynthesis */
    ac_am = assim(ci_am, gamma_star_am, vcmax_am, Km_am);
    ac_pm = assim(ci_pm, gamma_star_pm, vcmax_pm, Km_pm);

    /* Light-limited rate of photosynthesis allowed by RuBP regeneration */
    aj_am = assim(ci_am, gamma_star_am, jmax_am/4.0, 2.0*gamma_star_am);
    aj_pm = assim(ci_pm, gamma_star_pm, jmax_pm/4.0, 2.0*gamma_star_pm);

    /* light-saturated photosynthesis rate at the top of the canopy (gross) */
    asat_am = MIN(aj_am, ac_am);
    asat_pm = MIN(aj_pm, ac_pm);

    /* LUE (umol C umol-1 PAR) ; note conversion in epsilon */
    lue_am = epsilon(p, asat_am, m->par, alpha_am, daylen);
    lue_pm = epsilon(p, asat_pm, m->par, alpha_pm, daylen);

    /* use average to simulate canopy photosynthesis */
    lue_avg = (lue_am + lue_pm) / 2.0;

    /* absorbed photosynthetically active radiation (umol m-2 s-1) */
    if (float_eq(s->lai, 0.0))
        f->apar = 0.0;
    else
        f->apar = m->par * s->fipar;

    /* convert umol m-2 d-1 -> gC m-2 d-1 */
    conv = UMOL_TO_MOL * MOL_C_TO_GRAMS_C;
    f->gpp_gCm2 = f->apar * lue_avg * conv;
    f->gpp_am = (f->apar / 2.0) * lue_am * conv;
    f->gpp_pm = (f->apar / 2.0) * lue_pm * conv;

    /* g C m-2 to tonnes hectare-1 day-1 */
    f->gpp = f->gpp_gCm2 * G_AS_TONNES / M2_AS_HA;

    /* save apar in MJ m-2 d-1 */
    f->apar *= UMOL_2_JOL * J_TO_MJ;

    return;
}

double calculate_top_of_canopy_n(params *p, state *s, double ncontent)  {

    /*
    Calculate the canopy N at the top of the canopy (g N m-2), N0.
    Assuming an exponentially decreasing N distribution within the canopy:

    Returns:
    -------
    N0 : float (g N m-2)
        Top of the canopy N

    References:
    -----------
    * Chen et al 93, Oecologia, 93,63-69.

    */
    double N0;

    if (s->lai > 0.0) {
        /* calculation for canopy N content at the top of the canopy */
        N0 = ncontent * p->kext / (1.0 - exp(-p->kext * s->lai));
    } else {
        N0 = 0.0;
    }

    return (N0);
}

double calculate_co2_compensation_point(params *p, double Tk, double mt) {
    /*
        CO2 compensation point in the absence of mitochondrial respiration
        Rate of photosynthesis matches the rate of respiration and the net CO2
        assimilation is zero.

        Parameters:
        ----------
        Tk : float
            air temperature (Kelvin)

        Returns:
        -------
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
    */
    return (arrh(mt, p->gamstar25, p->eag, Tk));
}

double arrh(double mt, double k25, double Ea, double Tk) {
    /*
        Temperature dependence of kinetic parameters is described by an
        Arrhenius function

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.
    */
    return (k25 * exp((Ea * (Tk - mt)) / (mt * RGAS * Tk)));
}


double peaked_arrh(double mt, double k25, double Ea, double Tk, double deltaS,
                   double Hd) {
    /*
        Temperature dependancy approximated by peaked Arrhenius eqn,
        accounting for the rate of inhibition at higher temperatures.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [J mol-1]
        Tk : float
            leaf temperature [deg K]
        deltaS : float
            entropy factor [J mol-1 K-1)
        Hd : float
            describes rate of decrease about the optimum temp [J mol-1]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.

    */
    double arg1, arg2, arg3;

    arg1 = arrh(mt, k25, Ea, Tk);
    arg2 = 1.0 + exp((mt * deltaS - Hd) / (mt * RGAS));
    arg3 = 1.0 + exp((Tk * deltaS - Hd) / (Tk * RGAS));


    return (arg1 * arg2 / arg3);
}

double calculate_michaelis_menten_parameter(params *p, double Tk, double mt) {
    /*
        Effective Michaelis-Menten coefficent of Rubisco activity

        Parameters:
        ----------
        Tk : float
            air temperature (Kelvin)

        Returns:
        -------
        Km : float
            Effective Michaelis-Menten constant for Rubisco catalytic activity

        References:
        -----------
        Rubisco kinetic parameter values are from:
        * Bernacchi et al. (2001) PCE, 24, 253-259.
        * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    */

    double Kc, Ko;

    /* Michaelis-Menten coefficents for carboxylation by Rubisco */
    Kc = arrh(mt, p->kc25, p->eac, Tk);

    /* Michaelis-Menten coefficents for oxygenation by Rubisco */
    Ko = arrh(mt, p->ko25, p->eao, Tk);

    /* return effective Michaelis-Menten coefficient for CO2 */
    return ( Kc * (1.0 + p->oi / Ko) ) ;

}
void calculate_jmax_and_vcmax(control *c, params *p, state *s, double Tk,
                              double N0, double *jmax, double *vcmax,
                              double mt) {
    /*
        Calculate the maximum RuBP regeneration rate for light-saturated
        leaves at the top of the canopy (Jmax) and the maximum rate of
        rubisco-mediated carboxylation at the top of the canopy (Vcmax).

        Parameters:
        ----------
        Tk : float
            air temperature (Kelvin)
        N0 : float
            leaf N

        Returns:
        --------
        jmax : float (umol/m2/sec)
            the maximum rate of electron transport at 25 degC
        vcmax : float (umol/m2/sec)
            the maximum rate of electron transport at 25 degC
    */
    double jmax25, vcmax25;

    *vcmax = 0.0;
    *jmax = 0.0;

    if (c->modeljm == 0) {
        *jmax = p->jmax;
        *vcmax = p->vcmax;
    } else if (c->modeljm == 1) {
        /* the maximum rate of electron transport at 25 degC */
        jmax25 = p->jmaxna * N0 + p->jmaxnb;

        /* this response is well-behaved for TLEAF < 0.0 */
        *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk,
                            p->delsj, p->edj);

        /* the maximum rate of electron transport at 25 degC */
        vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
        *vcmax = vcmax_temperature(p, vcmax25, Tk - DEG_TO_KELVIN);
    } else if (c->modeljm == 2) {
        vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
        *vcmax = vcmax_temperature(p, vcmax25, Tk - DEG_TO_KELVIN);

        jmax25 = p->jv_slope * vcmax25 - p->jv_intercept;
        *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk, p->delsj,
                               p->edj);
    } else if (c->modeljm == 3) {
        /* the maximum rate of electron transport at 25 degC */
        jmax25 = p->jmax;

        /* this response is well-behaved for TLEAF < 0.0 */
        *jmax = peaked_arrh(mt, jmax25, p->eaj, Tk,
                            p->delsj, p->edj);

        /* the maximum rate of electron transport at 25 degC */
        vcmax25 = p->vcmax;
        *vcmax = vcmax_temperature(p, vcmax25, Tk - DEG_TO_KELVIN);

    }


    /* water stress acts through g1 (calculate_ci) and, optionally, on
       Jmax/Vcmax as a proxy for non-stomatal limitation */
    if (c->nonstomatal_limitation) {
        *jmax *= s->wtfac_root;
        *vcmax *= s->wtfac_root;
    }
    /*  Function allowing Jmax/Vcmax to be forced linearly to zero at low T */
    adj_for_low_temp(p, jmax, Tk);
    adj_for_low_temp(p, vcmax, Tk);

    return;

}

void adj_for_low_temp(params *p, double *param, double Tk) {
    /*
    Function allowing Jmax/Vcmax to be forced linearly to zero at low T

    Parameters:
    ----------
    Tk : float
        air temperature (Kelvin)
    */
    double lower_bound = p->photo_tlow;
    double upper_bound = p->photo_thigh;
    double Tc;

    Tc = Tk - DEG_TO_KELVIN;

    if (Tc < lower_bound)
        *param = 0.0;
    else if (Tc < upper_bound)
        *param *= (Tc - lower_bound) / (upper_bound - lower_bound);

    return;
}

double calculate_ci(control *c, params *p, state *s, double vpd, double Ca) {
    /*
    Calculate the intercellular (Ci) concentration

    Formed by substituting gs = g0 + 1.6 * (1 + (g1/sqrt(D))) * A/Ca into
    A = gs / 1.6 * (Ca - Ci) and assuming intercept (g0) = 0.

    Parameters:
    ----------
    vpd : float
        vapour pressure deficit [Pa]
    Ca : float
        ambient co2 concentration

    Returns:
    -------
    ci:ca : float
        ratio of intercellular to atmospheric CO2 concentration

    References:
    -----------
    * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    */

    double g1w, cica, ci=0.0;

    if (c->gs_model == MEDLYN) {
        g1w = p->g1 * s->wtfac_root;
        cica = g1w / (g1w + sqrt(vpd * PA_2_KPA));
        ci = cica * Ca;
    } else {
        prog_error("Only Belindas gs model is implemented", __LINE__);
    }

    return (ci);
}

double calculate_quantum_efficiency(params *p, double ci, double gamma_star) {
    /*

    Quantum efficiency for AM/PM periods replacing Sands 1996
    temperature dependancy function with eqn. from Medlyn, 2000 which is
    based on McMurtrie and Wang 1993.

    Parameters:
    ----------
    ci : float
        intercellular CO2 concentration.
    gamma_star : float [am/pm]
        CO2 compensation point in the abscence of mitochondrial respiration

    Returns:
    -------
    alpha : float
        Quantum efficiency

    References:
    -----------
    * Medlyn et al. (2000) Can. J. For. Res, 30, 873-888
    * McMurtrie and Wang (1993) PCE, 16, 1-13.

    */
    return (assim(ci, gamma_star, p->alpha_j/4.0, 2.0*gamma_star));
}

double assim(double ci, double gamma_star, double a1, double a2) {
    /*
    Morning and afternoon calcultion of photosynthesis with the
    limitation defined by the variables passed as a1 and a2, i.e. if we
    are calculating vcmax or jmax limited.

    Parameters:
    ----------
    ci : float
        intercellular CO2 concentration.
    gamma_star : float
        CO2 compensation point in the abscence of mitochondrial respiration
    a1 : float
        variable depends on whether the calculation is light or rubisco
        limited.
    a2 : float
        variable depends on whether the calculation is light or rubisco
        limited.

    Returns:
    -------
    assimilation_rate : float
        assimilation rate assuming either light or rubisco limitation.
    */
    if (ci < gamma_star)
        return (0.0);
    else
        return (a1 * (ci - gamma_star) / (a2 + ci));

}

double epsilon(params *p, double asat, double par, double alpha,
               double daylen) {
    /*
    Canopy scale LUE using method from Sands 1995, 1996.

    Sands derived daily canopy LUE from Asat by modelling the light response
    of photosysnthesis as a non-rectangular hyperbola with a curvature
    (theta) and a quantum efficiency (alpha).

    Assumptions of the approach are:
     - horizontally uniform canopy
     - PAR varies sinusoidally during daylight hours
     - extinction coefficient is constant all day
     - Asat and incident radiation decline through the canopy following
       Beer's Law.
     - leaf transmission is assumed to be zero.

    * Numerical integration of "g" is simplified to 6 intervals.

    Parameters:
    ----------
    asat : float
        Light-saturated photosynthetic rate at the top of the canopy
    par : float
        photosyntetically active radiation (umol m-2 d-1)
    theta : float
        curvature of photosynthetic light response curve
    alpha : float
        quantum yield of photosynthesis (mol mol-1)

    Returns:
    -------
    lue : float
        integrated light use efficiency over the canopy (umol C umol-1 PAR)

    Notes:
    ------
    NB. I've removed solar irradiance to PAR conversion. Sands had
    gamma = 2000000 to convert from SW radiation in MJ m-2 day-1 to
    umol PAR on the basis that 1 MJ m-2 = 2.08 mol m-2 & mol to umol = 1E6.
    We are passing PAR in umol m-2 d-1, thus avoiding the above.

    References:
    -----------
    See assumptions above...
    * Sands, P. J. (1995) Australian Journal of Plant Physiology,
      22, 601-14.

    */
    double delta, q, integral_g, sinx, arg1, arg2, arg3, lue, h;
    int i;

    /* subintervals scalar, i.e. 6 intervals */
    delta = 0.16666666667;

    /* number of seconds of daylight */
    h = daylen * SECS_IN_HOUR;

    if (asat > 0.0) {
        /* normalised daily irradiance */
        q = M_PI * p->kext * alpha * par / (2.0 * h * asat);
        integral_g = 0.0;
        for (i = 1; i < 13; i+=2) {
            sinx = sin(M_PI * i / 24.);
            arg1 = sinx;
            arg2 = 1.0 + q * sinx;
            arg3 = (sqrt(pow((1.0 + q * sinx), 2) - 4.0 * p->theta * q * sinx));
            integral_g += arg1 / (arg2 + arg3);
        }
        integral_g *= delta;
        lue = alpha * integral_g * M_PI;
    } else {
        lue = 0.0;
    }

    return (lue);
}


void mate_C4_photosynthesis(control *c, fluxes *f, met *m, params *p, state *s,
                            double daylen, double ncontent) {
    /*

    MATE simulates big-leaf C3 photosynthesis (GPP) based on Sands (1995),
    accounting for diurnal variations in irradiance and temp (am [sunrise-noon],
    pm[noon to sunset]).


    References:
    -----------
    * Collatz, G, J., Ribas-Carbo, M. and Berry, J. A. (1992) Coupled
      Photosynthesis-Stomatal Conductance Model for Leaves of C4 plants.
      Aust. J. Plant Physiol., 19, 519-38.
    * von Caemmerer, S. (2000) Biochemical Models of Leaf Photosynthesis. Chp 4.
      Modelling C4 photosynthesis. CSIRO PUBLISHING, Australia. pg 91-122.
    * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.

    Temperature dependancies:
    * Massad, R-S., Tuzet, A. and Bethenod, O. (2007) The effect of temperature
      on C4-type leaf photosynthesis parameters. Plant, Cell and Environment,
      30, 1191-1204.

    Intrinsic Quantum efficiency (mol mol-1), no Ci or temp dependancey
    in c4 plants see:
    * Ehleringer, J. R., 1978, Oecologia, 31, 255-267 or Collatz 1998.
    * Value taken from Table 1, Collatz et al.1998 Oecologia, 114, 441-454.

    */
    double N0, vcmax_am, vcmax_pm,
           ci_am, ci_pm, asat_am, asat_pm, lue_am, lue_pm, lue_avg,
           vcmax25_am, vcmax25_pm, par_per_sec, M_am, M_pm, A_am, A_pm,
           Rd_am, Rd_pm, conv;

    double mt = p->measurement_temp + DEG_TO_KELVIN;

    /* curvature parameter, transition between light-limited and
       carboxylation limited flux. Collatz table 2 */
    double beta1 = 0.83;

    /* curvature parameter, co-limitaiton between flux determined by
       Rubisco and light and CO2 limited flux. Collatz table 2 */
    double beta2 = 0.93;

    /* initial slope of photosynthetic CO2 response (mol m-2 s-1),
       Collatz table 2 */
    double kslope = 0.7;

    /* Calculate mate params & account for temperature dependencies */
    N0 = calculate_top_of_canopy_n(p, s, ncontent);

    ci_am = calculate_ci(c, p, s, m->vpd_am, m->Ca);
    ci_pm = calculate_ci(c, p, s, m->vpd_pm, m->Ca);

    /* Temp dependancies from Massad et al. 2007 */
    calculate_vcmax_parameter(c, p, s, m->Tk_am, N0, &vcmax_am, &vcmax25_am,
                              mt);
    calculate_vcmax_parameter(c, p, s, m->Tk_pm, N0, &vcmax_pm, &vcmax25_pm,
                              mt);

    /* Covert solar irradiance to PAR (umol PAR MJ-1) */
    conv = MJ_TO_J * J_2_UMOL;
    m->par *= conv;
    par_per_sec = m->par / (60.0 * 60.0 * daylen);

    /* Rubisco and light-limited capacity (Appendix, 2B) */
    M_am = quadratic(beta1, -(vcmax_am + p->alpha_c4 * par_per_sec),
                    (vcmax_am * p->alpha_c4 * par_per_sec));
    M_pm = quadratic(beta1, -(vcmax_pm + p->alpha_c4 * par_per_sec),
                    (vcmax_pm * p->alpha_c4 * par_per_sec));

    /* The limitation of the overall rate by M and CO2 limited flux */
    A_am = quadratic(beta2, -(M_am + kslope * ci_am),
                    (M_am * kslope * ci_am));
    A_pm = quadratic(beta2, -(M_pm + kslope * ci_pm),
                    (M_pm * kslope * ci_pm));

    /* These respiration terms are just for assimilation calculations,
       autotrophic respiration is stil assumed to be half of GPP */
    Rd_am = calc_respiration(m->Tk_am, vcmax25_am);
    Rd_pm = calc_respiration(m->Tk_pm, vcmax25_pm);

    /* Net (saturated) photosynthetic rate, not sure if this makes sense. */
    asat_am = A_am - Rd_am;
    asat_pm = A_pm - Rd_pm;

    /* LUE (umol C umol-1 PAR) */
    lue_am = epsilon(p, asat_am, m->par, p->alpha_c4, daylen);
    lue_pm = epsilon(p, asat_pm, m->par, p->alpha_c4, daylen);

    /* use average to simulate canopy photosynthesis */
    lue_avg = (lue_am + lue_pm) / 2.0;

    /* absorbed photosynthetically active radiation (umol m-2 s-1) */
    if (float_eq(s->lai, 0.0))
        f->apar = 0.0;
    else
        f->apar = m->par * s->fipar;

    /* convert umol m-2 d-1 -> gC m-2 d-1 */
    conv = UMOL_TO_MOL * MOL_C_TO_GRAMS_C;
    f->gpp_gCm2 = f->apar * lue_avg * conv;
    f->gpp_am = (f->apar / 2.0) * lue_am * conv;
    f->gpp_pm = (f->apar / 2.0) * lue_pm * conv;

    /* g C m-2 to tonnes hectare-1 day-1 */
    f->gpp = f->gpp_gCm2 * G_AS_TONNES / M2_AS_HA;

    /* save apar in MJ m-2 d-1 */
    f->apar *= UMOL_2_JOL * J_TO_MJ;

    return;
}


void calculate_vcmax_parameter(control *c, params *p, state *s, double Tk,
                               double N0,
                               double *vcmax, double *vcmax25, double mt) {
    /* Calculate the maximum rate of rubisco-mediated carboxylation at the
    top of the canopy

    References:
    -----------
    * Massad, R-S., Tuzet, A. and Bethenod, O. (2007) The effect of
      temperature on C4-type leaf photosynthesis parameters. Plant, Cell and
      Environment, 30, 1191-1204.

    # http://www.cesm.ucar.edu/models/cesm1.0/clm/CLM4_Tech_Note.pdf
    # Table 8.2 has PFT values...

    Parameters:
    ----------
    Tk : float
        air temperature (kelvin)
    N0 : float
        leaf N

    Returns:
    -------
    vcmax : float, list [am, pm]
        maximum rate of Rubisco activity
    */
    /* Massad et al. 2007 */
    /*Ea = self.params.eav
    Hd = self.params.edv
    delS = self.params.delsv */
    double Ea = 67294.0;
    double Hd = 144568.0;
    double delS = 472.0;

    /* the maximum rate of electron transport at 25 degC */
    *vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
    *vcmax = peaked_arrh(mt, *(vcmax25), Ea, Tk, delS, Hd);

    /* water stress acts through g1 (calculate_ci) and, optionally, on
       Vcmax as a proxy for non-stomatal limitation */
    if (c->nonstomatal_limitation) {
        *vcmax *= s->wtfac_root;
    }

    /* Function allowing Jmax/Vcmax to be forced linearly to zero at low T */
    adj_for_low_temp(p, vcmax, Tk);

    return;
}

double calc_respiration(double Tk, double vcmax25) {
    /*
    Mitochondrial respiration may occur in the mesophyll as well as in the
    bundle sheath. As rubisco may more readily refix CO2 released in the
    bundle sheath, Rd is described by its mesophyll and bundle-sheath
    components: Rd = Rm + Rs

    Parameters:
    ----------
    Tk : float
        air temperature (kelvin)
    vcmax25 : float, list

    Returns:
    -------
    Rd : float, list [am, pm]
        (respiration in the light) 'day' respiration (umol m-2 s-1)

    References:
    -----------
    Tjoelker et al (2001) GCB, 7, 223-230.
    */
    double Rd25, Tref = 25.0;

    /* ratio between respiration rate at one temperature and the respiration
       rate at a temperature 10 deg C lower */
    double Q10 = 2.0;

    /* scaling constant to Vcmax25, value = 0.015 after Collatz et al 1991.
       Agricultural and Forest Meteorology, 54, 107-136. But this if for C3
       Value used in JULES for C4 is 0.025, using that one, see Clark et al.
       2011, Geosci. Model Dev, 4, 701-722. */
    double fdr = 0.025;

    /* specific respiration at a reference temperature (25 deg C) */
    Rd25 = fdr * vcmax25;

    return (Rd25 * pow(Q10, (((Tk - DEG_TO_KELVIN) - Tref) / 10.0)));
}

double quadratic(double a, double b, double c) {
    /* minimilist quadratic solution

    Parameters:
    ----------
    a : float
        co-efficient
    b : float
        co-efficient
    c : float
        co-efficient

    Returns:
    -------
    root : float

    */
    double d, root;

    /* discriminant */
    d = pow(b, 2.0) - 4.0 * a * c;

    /* Negative quadratic equation */
    root = (-b - sqrt(d)) / (2.0 * a);

    return (root);
}
