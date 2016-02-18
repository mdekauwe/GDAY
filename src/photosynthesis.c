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
*   14.01.2016
*
* =========================================================================== */
#include "photosynthesis.h"


void photosynthesis_C3(control *c, params *p, state *s, double N0,
                       double tleaf, double par, double Cs, double dleaf,
                       double *anleaf, double *gsc, bool leaf) {
    /*
        Calculate photosynthesis following Farquhar & von Caemmerer, this is an
        implementation of the routinue in MAESTRA

        References:
        -----------
        * GD Farquhar, S Von Caemmerer (1982) Modelling of photosynthetic
          response to environmental conditions. Encyclopedia of plant
          physiology 12, 549-587.
        * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
        * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.
    */

    double gamma_star, km, jmax, vcmax, rd, J, Vj, gs_over_a, g0;
    double A, B, C, arg1, arg2, Ci, Ac, Aj, dleaf_kpa;
    double Rd0 = 0.92; /* Dark respiration rate make a paramater! */
    int    qudratic_error = FALSE, large_root;
    double g0_zero = 1E-09; /* numerical issues, don't use zero */

    /* Calculate photosynthetic parameters from leaf temperature. */
    gamma_star = calc_co2_compensation_point(p, tleaf);
    km = calculate_michaelis_menten(p, tleaf);
    calculate_jmaxt_vcmaxt(c, p, s, tleaf, N0, &jmax, &vcmax, leaf);

    /* leaf respiration in the light, Collatz et al. 1991 */
    rd = 0.015 * vcmax;
    /*rd = calc_leaf_day_respiration(tleaf, Rd0); */

    /* actual electron transport rate */
    qudratic_error = FALSE;
    large_root = FALSE;
    J = quad(p->theta, -(p->alpha_j * par + jmax),
             p->alpha_j * par * jmax, large_root, &qudratic_error);
    /* RuBP regeneration rate */
    Vj = J / 4.0;

    /* Deal with extreme cases */
    if (jmax <= 0.0 || vcmax <= 0.0 || isnan(J)) {
        *anleaf = -rd;
        *gsc = g0_zero;
    } else {
        /* Hardwiring this for Medlyn gs model for the moment, till I figure
        out the best structure */

        /* For the medlyn model this is already in conductance to CO2, so the
           1.6 from the corrigendum to Medlyn et al 2011 is missing here */
        dleaf_kpa = dleaf * PA_2_KPA;
        if (dleaf_kpa < 0.05) {
            dleaf_kpa = 0.05;
        }
        gs_over_a = (1.0 + (p->g1 * s->wtfac_root) / sqrt(dleaf_kpa)) / Cs;
        g0 = g0_zero;

        /* Solution when Rubisco activity is limiting */
        A = g0 + gs_over_a * (vcmax - rd);
        B = ( (1.0 - Cs * gs_over_a) * (vcmax - rd) + g0 * (km - Cs) -
               gs_over_a * (vcmax * gamma_star + km * rd) );
        C = ( -(1.0 - Cs * gs_over_a) * (vcmax * gamma_star + km * rd) -
               g0 * km * Cs );

        /* intercellular CO2 concentration */
        qudratic_error = FALSE;
        large_root = TRUE;
        Ci = quad(A, B, C, large_root, &qudratic_error);

        if (qudratic_error || Ci <= 0.0 || Ci > Cs) {
            Ac = 0.0;
        } else {
            Ac = vcmax * (Ci - gamma_star) / (Ci + km);
        }

        /* Solution when electron transport rate is limiting */
        A = g0 + gs_over_a * (Vj - rd);
        B = ( (1. - Cs * gs_over_a) * (Vj - rd) + g0 *
              (2. * gamma_star - Cs) - gs_over_a *
              (Vj * gamma_star + 2.0 * gamma_star * rd) );
        C = ( -(1.0 - Cs * gs_over_a) * gamma_star * (Vj + 2.0 * rd) -
               g0 * 2.0 * gamma_star * Cs );

        /* intercellular CO2 concentration */
        qudratic_error = FALSE;
        large_root = TRUE;
        Ci = quad(A, B, C, large_root, &qudratic_error);

        Aj = Vj * (Ci - gamma_star) / (Ci + 2.0 * gamma_star);
        /* Below light compensation point? */
        if (Aj - rd < 1E-6) {
            Ci = Cs;
            Aj = Vj * (Ci - gamma_star) / (Ci + 2.0 * gamma_star);
        }
        *anleaf = MIN(Ac, Aj) - rd;
        *gsc = MAX(g0, g0 + gs_over_a * *anleaf);
    }
    return;
}

double calc_co2_compensation_point(params *p, double tleaf) {
    /*
        CO2 compensation point in the absence of non-photorespiratory
        respiration.

        Parameters:
        ----------
        tleaf : float
            leaf temperature (deg C)

        Returns:
        -------
        gamma_star : float
            CO2 compensation point in the abscence of mitochondrial respiration
            (umol mol-1)

        References:
        -----------
        * Bernacchi et al 2001 PCE 24: 253-260
    */
    return (arrhenius(p->gamstar25, p->eag, tleaf, p->measurement_temp));
}

double calculate_michaelis_menten(params *p, double tleaf) {
    /*
        Effective Michaelis-Menten coefficent of Rubisco activity

        Parameters:
        ----------
        tleaf : float
            leaf temperature (deg C)

        Returns:
        -------
        Km : float
            Effective Michaelis-Menten constant for Rubisco catalytic activity
            (umol mol-1)

        References:
        -----------
        Rubisco kinetic parameter values are from:
        * Bernacchi et al. (2001) PCE, 24, 253-259.
        * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    */
    double Kc, Ko, Km;

    /* Michaelis-Menten coefficents for carboxylation by Rubisco */
    Kc = arrhenius(p->kc25, p->eac, tleaf, p->measurement_temp);

    /* Michaelis-Menten coefficents for oxygenation by Rubisco */
    Ko = arrhenius(p->ko25, p->eao, tleaf, p->measurement_temp);

    /* return effective Michaelis-Menten coefficient for CO2 */
    Km = Kc * (1.0 + p->oi / Ko);

    return (Km);
}

void calculate_jmaxt_vcmaxt(control *c, params *p, state *s, double tleaf,
                            double N0, double *jmax, double *vcmax,
                            bool leaf) {
    /*
        Calculate the potential electron transport rate (Jmax) and the
        maximum Rubisco activity (Vcmax) at the leaf temperature.

        For Jmax -> peaked arrhenius is well behaved for tleaf < 0.0

        Parameters:
        ----------
        tleaf : float
            air temperature (deg C)
        N0 : float
            leaf N at the top of the canopy
        jmax : float
            the potential electron transport rate at the leaf temperature
            (umol m-2 s-1)
        vcmax : float
            the maximum Rubisco activity at the leaf temperature (umol m-2 s-1)
    */
    double jmax25, vcmax25, a, b;
    double lower_bound = 0.0;
    double upper_bound = 10.0;
    double tref = p->measurement_temp;
    *vcmax = 0.0;
    *jmax = 0.0;

    if (c->modeljm == 0) {
        *jmax = p->jmax;
        *vcmax = p->vcmax;
    } else if (c->modeljm == 1) {

        if (leaf == SUNLIT) {
            jmax25 = integrate_sunlit_frac(p->jmaxna, p->jmaxnb, N0, s->lai);
            vcmax25 = integrate_sunlit_frac(p->vcmaxna, p->vcmaxnb, N0, s->lai);
        } else {
            jmax25 = integrate_shaded_frac(p->jmaxna, p->jmaxnb, N0, s->lai);
            vcmax25 = integrate_shaded_frac(p->vcmaxna, p->vcmaxnb, N0,s->lai);
        }
        *jmax = peaked_arrhenius(jmax25, p->eaj, tleaf, tref, p->delsj, p->edj);
        *vcmax = arrhenius(vcmax25, p->eav, tleaf, tref);

        printf("%d %lf %lf\n", leaf, *vcmax, *jmax);
        /*
        jmax25 = p->jmaxna * N0 + p->jmaxnb;
        *jmax = peaked_arrhenius(jmax25, p->eaj, tleaf, tref, p->delsj, p->edj);
        vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
        *vcmax = arrhenius(vcmax25, p->eav, tleaf, tref);*/
    } else if (c->modeljm == 2) {
        vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
        *vcmax = arrhenius(vcmax25, p->eav, tleaf, tref);
        jmax25 = p->jv_slope * vcmax25 - p->jv_intercept;
        *jmax = peaked_arrhenius(jmax25, p->eaj, tleaf, tref, p->delsj, p->edj);
    } else if (c->modeljm == 3) {
        jmax25 = p->jmax;
        *jmax = peaked_arrhenius(jmax25, p->eaj, tleaf, tref, p->delsj, p->edj);
        vcmax25 = p->vcmax;
        *vcmax = arrhenius(vcmax25, p->eav, tleaf, tref);
    }

    /* reduce photosynthetic capacity with moisture stress */
    *jmax *= s->wtfac_root;
    *vcmax *= s->wtfac_root;

    /* Jmax/Vcmax forced linearly to zero at low T */
    if (tleaf < lower_bound) {
        *jmax = 0.0;
        *vcmax = 0.0;
    } else if (tleaf < upper_bound) {
        *jmax *= (tleaf - lower_bound) / (upper_bound - lower_bound);
        *vcmax *= (tleaf - lower_bound) / (upper_bound - lower_bound);
    }

    return;
}

double integrate_sunlit_frac(double a, double b, double N0,
                                     double lai) {

    /*
        Integrate over the canopy depth to yeild bulk values of sunlit
        vcmax/jmax

        References:
        ----------
        * Dai et al. (2004) Journal of Climate, 17, 2281-2299., eqn 37a,b
    */
    double kb = 0.5;
    double kn = 0.3;
    double sun, arg1, arg2;

    arg1 = a / kb * (1.0 - exp(-kb * lai));
    arg2 = b * N0 / (kb + kn) * (1.0 - exp(-(kn + kb) * lai));
    sun = arg1 + arg2;

    return (sun);
}

double integrate_shaded_frac(double a, double b, double N0,
                                     double lai) {

    /*
        Integrate over the canopy depth to yeild bulk values of shaded
        vcmax/jmax

        References:
        ----------
        * Dai et al. (2004) Journal of Climate, 17, 2281-2299., eqn 38a,b
    */
    double kb = 0.5;
    double kn = 0.3;
    double shade, arg1, arg2, arg3, arg4;

    arg1 = a * lai;
    arg2 = (a / kb) * (exp(-kb * lai) - 1.0);
    arg3 = b * N0 / kn * (exp(-kb * lai) - 1.0);
    arg4 = b * N0 / (kn + kb) * (exp(-(kn + kb) * lai) - 1.0);
    shade = arg1 + arg2 - arg2 + arg3;

    return (shade);
}



double calc_leaf_day_respiration(double tleaf, double Rd0) {
    /* Calculate leaf respiration in the light using a Q10 (exponential)
    formulation

    Parameters:
    ----------
    Rd0 : float
        Estimate of respiration rate at the reference temperature 25 deg C
        or or 298 K
    Tleaf : float
        leaf temperature

    Returns:
    -------
    Rd : float
        leaf respiration in the light

    */
    double Rd;
    /* tbelow is the minimum temperature at which respiration occurs;
       below this, Rd(T) is set to zero. */
    double tbelow = -100.0;

    /* Determines the respiration in the light, relative to
       respiration in the dark. For example, if DAYRESP = 0.7, respiration
       in the light is 30% lower than in the dark. */
    double day_resp = 1.0;  /* no effect of light */

    /* Default temp at which maint resp specified */
    double rtemp = 0.0;

    /* exponential coefficient of the temperature response of foliage
       respiration */
    double q10f = 0.067;

    /* Temperature sensitivity of Rd. */
    double Q10 = 2.0;

    if (tleaf >= tbelow)
        Rd = Rd0 * day_resp * exp(q10f * (tleaf - rtemp));
    else
        Rd = 0.0;

    return (Rd);
}

double arrhenius(double k25, double Ea, double T, double Tref) {
    /*
        Temperature dependence of kinetic parameters is described by an
        Arrhenius function

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [J mol-1]
        T : float
            temperature [deg C]
        Tref : float
            measurement temperature [deg C]

        Returns:
        -------
        kt : float
            temperature dependence on parameter

        References:
        -----------
        * Medlyn et al. 2002, PCE, 25, 1167-1179.
    */
    double Tk, TrefK;
    Tk = T + DEG_TO_KELVIN;
    TrefK = Tref + DEG_TO_KELVIN;

    return (k25 * exp(Ea * (T - Tref) / (RGAS * Tk * TrefK)));
}

double peaked_arrhenius(double k25, double Ea, double T, double Tref,
                        double deltaS, double Hd) {
    /*
        Temperature dependancy approximated by peaked Arrhenius eqn,
        accounting for the rate of inhibition at higher temperatures.

        Parameters:
        ----------
        k25 : float
            rate parameter value at 25 degC
        Ea : float
            activation energy for the parameter [J mol-1]
        T : float
            temperature [deg C]
        Tref : float
            measurement temperature [deg C]
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
    double Tk, TrefK;
    Tk = T + DEG_TO_KELVIN;
    TrefK = Tref + DEG_TO_KELVIN;

    arg1 = arrhenius(k25, Ea, T, Tref);
    arg2 = 1.0 + exp((deltaS * TrefK - Hd) / (RGAS * TrefK));
    arg3 = 1.0 + exp((deltaS * Tk - Hd) / (RGAS * Tk));

    return (arg1 * arg2 / arg3);
}


double quad(double a, double b, double c, bool large, int *error) {
    /* quadratic solution

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
    d = (b * b) - 4.0 * a * c;
    if (d < 0.0) {
        /*fprintf(stderr, "imaginary root found\n");
        exit(EXIT_FAILURE);*/
        *error = TRUE;
    }

    if (large) {
        if (a == 0.0 && b > 0.0) {
            root = -c / b;
        } else if (a == 0.0 && b == 0.0) {
            root = 0.0;
            if (c != 0.0) {
                /*fprintf(stderr, "Can't solve quadratic\n");
                exit(EXIT_FAILURE);*/
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
                /*fprintf(stderr, "Can't solve quadratic\n");
                exit(EXIT_FAILURE);*/
                *error = TRUE;
            }
        } else {
            root = (-b - sqrt(d)) / (2.0 * a);
        }
    }
    return (root);
}








void mate_C3_photosynthesis(control *c, fluxes *f, met *m, params *p, state *s,
                            int project_day, double daylen, double ncontent) {
    /*

    MATE simulates big-leaf C3 photosynthesis (GPP) based on Sands (1995),
    accounting for diurnal variations in irradiance and temp (am [sunrise-noon],
    pm[noon to sunset]).

    MATE is connected to G'DAY via LAI and leaf N content.

    Plant autotrophic respiration is calculated via carbon-use efficiency
    (CUE=NPP/GPP).

    References:
    -----------
    * Medlyn, B. E. et al (2011) Global Change Biology, 17, 2134-2144.
    * McMurtrie, R. E. et al. (2008) Functional Change Biology, 35, 521-34.
    * Sands, P. J. (1995) Australian Journal of Plant Physiology, 22, 601-14.

    Rubisco kinetic parameter values are from:
    * Bernacchi et al. (2001) PCE, 24, 253-259.
    * Medlyn et al. (2002) PCE, 25, 1167-1179, see pg. 1170.

    */
    double N0, Tk_am, Tk_pm, par, vpd_am, vpd_pm, ca, gamma_star_am,
           gamma_star_pm, Km_am, Km_pm, jmax_am, jmax_pm, vcmax_am, vcmax_pm,
           ci_am, ci_pm, alpha_am, alpha_pm, ac_am, ac_pm, aj_am, aj_pm,
           asat_am, asat_pm, lue_am, lue_pm, lue_avg, conv;
    double mt = p->measurement_temp + DEG_TO_KELVIN;

    get_met_stuff(m, project_day, &Tk_am, &Tk_pm, &par, &vpd_am,
                  &vpd_pm, &ca);

    /* Calculate mate params & account for temperature dependencies */
    N0 = calculate_top_of_canopy_n(p, s, ncontent);

    gamma_star_am = calculate_co2_compensation_point(p, Tk_am, mt);
    gamma_star_pm = calculate_co2_compensation_point(p, Tk_pm, mt);

    Km_am = calculate_michaelis_menten_parameter(p, Tk_am, mt);
    Km_pm = calculate_michaelis_menten_parameter(p, Tk_pm, mt);

    calculate_jmax_and_vcmax(c, p, s, Tk_am, N0, &jmax_am, &vcmax_am, mt);
    calculate_jmax_and_vcmax(c, p, s, Tk_pm, N0, &jmax_pm, &vcmax_pm, mt);

    ci_am = calculate_ci(c, p, s, vpd_am, ca);
    ci_pm = calculate_ci(c, p, s, vpd_pm, ca);

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

    /* LUE (umol C umol-1 PAR) */
    lue_am = epsilon(p, asat_am, par, daylen, alpha_am);
    lue_pm = epsilon(p, asat_pm, par, daylen, alpha_pm);

    /* use average to simulate canopy photosynthesis */
    lue_avg = (lue_am + lue_pm) / 2.0;

    /* absorbed photosynthetically active radiation (umol m-2 s-1) */
    if (float_eq(s->lai, 0.0))
        f->apar = 0.0;
    else
        f->apar = par * s->fipar;

    /* convert umol m-2 s-1 -> gC m-2 d-1 */
    conv = UMOL_TO_MOL * MOL_C_TO_GRAMS_C * (60.0 * 60.0 * daylen);
    f->gpp_gCm2 = f->apar * lue_avg * conv;
    f->gpp_am = (f->apar / 2.0) * lue_am * conv;
    f->gpp_pm = (f->apar / 2.0) * lue_pm * conv;

    /* g C m-2 to tonnes hectare-1 day-1 */
    f->gpp = f->gpp_gCm2 * G_AS_TONNES / M2_AS_HA;

    return;
}

void get_met_stuff(met *m, int project_day, double *Tk_am,
                   double *Tk_pm, double *par, double *vpd_am, double *vpd_pm,
                   double *ca) {
    /*
    Grab the days met data out of the structure and return day values.

    Parameters:
    ----------
    project_day : int
        project day.

    Returns:
    -------
    TK_am : float
        am air temperature [Kelvin]
    TK_am : float
        pm air temperature [Kelvin]
    vpd_am : float
        am vpd [kPa]
    vpd_pm : float
        pm vpd [kPa]
    par : float
        average daytime PAR [umol m-2 s-1]
    ca : float
        atmospheric co2 [umol mol-1]
    */

    *Tk_am = m->tam[project_day] + DEG_TO_KELVIN;
    *Tk_pm = m->tpm[project_day] + DEG_TO_KELVIN;
    *vpd_am = m->vpd_am[project_day];
    *vpd_pm = m->vpd_pm[project_day];
    *ca = m->co2[project_day];
    *par = m->par[project_day];

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
        *vcmax = arrh(mt, vcmax25, p->eav, Tk);
    } else if (c->modeljm == 2) {
        vcmax25 = p->vcmaxna * N0 + p->vcmaxnb;
        *vcmax = arrh(mt, vcmax25, p->eav, Tk);

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
        *vcmax = arrh(mt, vcmax25, p->eav, Tk);

    }


    /* reduce photosynthetic capacity with moisture stress */
    *jmax *= s->wtfac_root;
    *vcmax *= s->wtfac_root;
    /*  Function allowing Jmax/Vcmax to be forced linearly to zero at low T */
    adj_for_low_temp(*(&jmax), Tk);
    adj_for_low_temp(*(&vcmax), Tk);

    return;

}

void adj_for_low_temp(double *param, double Tk) {
    /*
    Function allowing Jmax/Vcmax to be forced linearly to zero at low T

    Parameters:
    ----------
    Tk : float
        air temperature (Kelvin)
    */
    double lower_bound = 0.0;
    double upper_bound = 10.0;
    double Tc;

    Tc = Tk - DEG_TO_KELVIN;

    if (Tc < lower_bound)
        *param = 0.0;
    else if (Tc < upper_bound)
        *param *= (Tc - lower_bound) / (upper_bound - lower_bound);

    return;
}

double calculate_ci(control *c, params *p, state *s, double vpd, double ca) {
    /*
    Calculate the intercellular (Ci) concentration

    Formed by substituting gs = g0 + 1.6 * (1 + (g1/sqrt(D))) * A/Ca into
    A = gs / 1.6 * (Ca - Ci) and assuming intercept (g0) = 0.

    Parameters:
    ----------
    vpd : float
        vapour pressure deficit
    ca : float
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
        cica = g1w / (g1w + sqrt(vpd));
        ci = cica * ca;
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

double epsilon(params *p, double asat, double par, double daylen,
               double alpha) {
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
    daylen : float
        length of day (hrs).
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
    NB. I've removed some of the unit conversions as they are unnecessary.
    Sands had gamma = 2000000 to convert from SW radiation in MJ m-2 day-1 to
    umol PAR on the basis that 1 MJ m-2 = 2.08 mol m-2 & mol to umol = 1E6.
    He then has a divide by h (where h is in seconds) to go from day to seconds
    of PAR. But if we just pass PAR in umol m-2 s-1 we don't need all of this

    References:
    -----------
    See assumptions above...
    * Sands, P. J. (1995) Australian Journal of Plant Physiology,
      22, 601-14.

    */
    double delta, q, integral_g, sinx, arg1, arg2, arg3, lue;
    int i;

    /* subintervals scaler, i.e. 6 intervals */
    delta = 0.16666666667;

    if (asat > 0.0) {
        /* normalised daily irradiance */
        q = M_PI * p->kext * alpha * par / (2.0 * asat);
        integral_g = 0.0;
        for (i = 1; i < 13; i+=2) {
            sinx = sin(M_PI * i / 24.);
            arg1 = sinx;
            arg2 = 1.0 + q * sinx;
            arg3 = (sqrt(pow((1.0 + q * sinx), 2) - 4.0 * p->theta * q * sinx));
            integral_g += arg1 / (arg2 + arg3) * delta;
        }
        lue = alpha * integral_g * M_PI;
    } else {
        lue = 0.0;
    }

    return (lue);
}




void mate_C4_photosynthesis(control *c, fluxes *f, met *m, params *p, state *s,
                            int project_day, double daylen, double ncontent) {
    /*

    MATE simulates big-leaf C3 photosynthesis (GPP) based on Sands (1995),
    accounting for diurnal variations in irradiance and temp (am [sunrise-noon],
    pm[noon to sunset]).

    MATE is connected to G'DAY via LAI and leaf N content.

    Plant autotrophic respiration is calculated via carbon-use efficiency
    (CUE=NPP/GPP).

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
    double N0, Tk_am, Tk_pm, par, vpd_am, vpd_pm, ca, vcmax_am, vcmax_pm,
           ci_am, ci_pm, asat_am, asat_pm, lue_am, lue_pm, lue_avg,
           vcmax25_am, vcmax25_pm, par_per_sec, M_am, M_pm, A_am, A_pm,
           Rd_am, Rd_pm, apar_half_day, SEC_2_HALF_DAY, conv;

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

    get_met_stuff(m, project_day, &Tk_am, &Tk_pm, &par, &vpd_am,
                  &vpd_pm, &ca);

    /* Calculate mate params & account for temperature dependencies */
    N0 = calculate_top_of_canopy_n(p, s, ncontent);

    ci_am = calculate_ci(c, p, s, vpd_am, ca);
    ci_pm = calculate_ci(c, p, s, vpd_pm, ca);

    /* Temp dependancies from Massad et al. 2007 */
    calculate_vcmax_parameter(p, s, Tk_am, N0, &vcmax_am, &vcmax25_am, mt);
    calculate_vcmax_parameter(p, s, Tk_pm, N0, &vcmax_pm, &vcmax25_pm, mt);

    /* Rubisco and light-limited capacity (Appendix, 2B) */
    par_per_sec = par / (60.0 * 60.0 * daylen);
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
    Rd_am = calc_respiration(Tk_am, vcmax25_am);
    Rd_pm = calc_respiration(Tk_pm, vcmax25_pm);

    /* Net (saturated) photosynthetic rate, not sure if this makes sense. */
    asat_am = A_am - Rd_am;
    asat_pm = A_pm - Rd_pm;

    /* LUE (umol C umol-1 PAR) */
    lue_am = epsilon(p, asat_am, par, daylen, p->alpha_c4);
    lue_pm = epsilon(p, asat_pm, par, daylen, p->alpha_c4);

    /* use average to simulate canopy photosynthesis */
    lue_avg = (lue_am + lue_pm) / 2.0;

    /* absorbed photosynthetically active radiation (umol m-2 s-1) */
    if (float_eq(s->lai, 0.0))
        f->apar = 0.0;
    else
        f->apar = par * s->fipar;

    /* convert umol m-2 d-1 -> gC m-2 d-1 */
    conv = UMOL_2_GRAMS_C * (60.0 * 60.0 * daylen);
    f->gpp_gCm2 = f->apar * lue_avg * conv;
    f->gpp_am = (f->apar / 2.0) * lue_am * conv;
    f->gpp_pm = (f->apar / 2.0) * lue_pm * conv;
    f->npp_gCm2 = f->gpp_gCm2 * p->cue;

    /* g C m-2 to tonnes hectare-1 day-1 */
    f->gpp = f->gpp_gCm2 * GRAM_C_2_TONNES_HA;
    f->npp = f->npp_gCm2 * GRAM_C_2_TONNES_HA;

    return;
}


void calculate_vcmax_parameter(params *p, state *s, double Tk, double N0,
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

    /* reduce photosynthetic capacity with moisture stress */
    *vcmax *= s->wtfac_root;

    /* Function allowing Jmax/Vcmax to be forced linearly to zero at low T */
    adj_for_low_temp(*(&vcmax), Tk);

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














/*


        # Reducing assimilation if we encounter frost. Frost is assumed to
        # impact on the maximum photosynthetic capacity and alpha_j
        # So there is only an indirect effect on LAI, this could be changed...
        if c->frost:
            Tmax = m->'tmax'][day]
            Tmin = m->'tmin'][day]

            Thard = self.calc_frost_hardiness(daylen, Tmin, Tmax)
            (total_alpha_limf,
            total_amax_limf) = self.calc_frost_impact_factors(Thard, Tmin, Tmax)
            alpha_am *= total_alpha_limf
            alpha_pm *= total_alpha_limf





]


    def calc_frost_hardiness(self, daylength, Tmin, Tmax):
        """ Capacity of plants to survive frost defined by a hardiness
        paramater, Thard.


        Parameters:
        ----------


        Returns:
        -------
        Thard : float (deg C)
            Nightly minimum temperature causing a 50% reduction in Amax in
            previously undamaged leaves.

        References:
        -----------
        * King and Ball, 1998, Aust. J. Plant Physiol., 25, 27-37.
        """
        beta = 1.0 # degC/h

        # Average night-time temperature
        Tnight = Tmin + 0.25 * (Tmax - Tmin)

        # equinox daylength
        Teq = 12.0

        # Stationary level of hardiness
        Tstat = (p->frost_a + p->frost_b *
                 (Tnight + beta * (daylength - Teq)))

        # previous days Thard
        if p->thardp is None:
            p->thardp = Tstat

        # Frost hardiness parameter
        Thard = (p->thardp + p->frost_c *
                 (Tstat - p->thardp))

        if Thard < -12.0:
            Thard = -12.0
        elif Thard > -3.0:
            Thard = -3.0

        # set previous value to todays value
        p->thardp = Thard

        return (Thard)

    def calc_frost_impact_factors(self, Thard, Tmin, Tmax):
        """ Calculate multiplicative frost impact factors, 0=lethal frost; 1=no
        damage from the previous night


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
        total_alpha_limf : float [0-1]
            limitation on alpha
        total_amax_limf : float [0-1]
            limitation on Amax

        References:
        -----------
        * King and Ball, 1998, Aust. J. Plant Physiol., 25, 27-37.
        """

        # Temperature range between 0 and 100% photosynthetic damage from low
        # temp following Battaglia et a. 2004.
        Trange = 5.0
        #Trange = Tmax - Tmin

        # Factor accounting for the previous nights frost on Amax
        if Tmin > Thard + 0.5 * Trange:
            f_A = 1.0
        elif Thard + 0.5 * Trange > Tmin and Tmin > Thard - 0.5 * Trange:
            f_A = 0.5 * (1.0 + sin(pi * (Tmin - Thard) / Trange))
        elif Tmin <= Thard - 0.5 * Trange:
            f_A = 0.0

        # Factor accounting for effect on initial slope of the light response
        # curve (alpha)
        d = 1.5 + 0.5 * exp(-2.0 * p->kext * s->lai)
        if f_A >= 0.5:
            f_alpha = f_A**d
        elif f_A < 0.5:
            f_alpha = f_A / d

        #
        ## Short term effects, complete recovery ~ 5 days
        #
        if p->fcap < 0.8:
            fcA = f_A * (p->fcap + 0.2)
        elif p->fcap >= 0.8:
            fcA = f_A

        if p->fc_alpha_p < 0.8:
            fc_alpha = f_alpha * (p->fc_alpha_p + 0.2)
        elif p->fc_alpha_p >= 0.8:
            fc_alpha = f_alpha

        # set previous days values to todays value
        p->fcap = fcA
        p->fc_alpha_p = fc_alpha

        #
        ## Long term cumulative frost impact factor
        #
        if f_alpha < 1.0:
            f_long = f_alpha**p->frost_p * p->f_long_gp
        elif f_alpha == 1.0:
            f_long = 0.01 + 0.99 * p->f_long_gp

        # set previous value to todays value
        p->f_long_gp = f_long

        # Combined factor
        total_alpha_limf = max(min(1.0, f_long * fc_alpha), 0.0)
        total_amax_limf = max(min(1.0, fcA * f_long), 0.0)

        return (total_alpha_limf, total_amax_limf)
*/
