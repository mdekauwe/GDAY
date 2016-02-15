#include "water_balance.h"


void calculate_water_balance(control *c, fluxes *f, met *m, params *p,
                             state *s, int day_idx, int daylen,
                             double trans_leaf, double omega_leaf) {
    /*

    Calculate the water balance (including all water fluxes).

    Parameters:
    ----------
    control : structure
        control structure
    fluxes : structure
        fluxes structure
    met : structure
        meteorological drivers structure
    params : structure
        parameters structure
    day : int
        project day. (Dummy argument, only passed for daily model)
    daylen : double
        length of day in hours. (Dummy argument, only passed for daily model)
    trans_leaf : double
        leaf transpiration (Dummy argument, only passed for sub-daily model)
    omega_leaf : double
        decoupling coefficient (Dummy argument, only passed for sub-daily model)

    */
    double soil_evap, et, interception, runoff, rain, press, sw_rad, conv,
           tair, transpiration, net_rad, SEC_2_HALF_DAY, HALF_DAY_2_SEC,
           transpiration_am, transpiration_pm, gs_am, gs_pm, LE_am,
           LE_pm, ga_am, ga_pm, net_rad_day, net_rad_am, net_rad_pm, trans_am,
           omega_am, gs_mol_m2_hfday_am, ga_mol_m2_hfday_am, tair_am, tair_pm,
           tair_day, sw_rad_am, sw_rad_pm, sw_rad_day, vpd_am, vpd_pm, vpd_day,
           wind_am, wind_pm, wind_day, ca, gpp_am, gpp_pm, trans_pm,
           omega_pm, gs_mol_m2_hfday_pm, ga_mol_m2_hfday_pm;

    SEC_2_HALF_DAY =  60.0 * 60.0 * (daylen / 2.0);
    HALF_DAY_2_SEC = 1.0 / SEC_2_HALF_DAY;

    /* unpack met forcing */
    if (c->sub_daily) {
        rain = m->rain[c->hrly_idx];
        press = m->press[c->hrly_idx] * KPA_2_PA;
        tair = m->tair[c->hrly_idx];
        sw_rad = m->par[c->hrly_idx] * PAR_2_SW; /* W m-2 */
    } else {
        ca = m->co2[day_idx];
        tair = m->tair[day_idx];
        tair_am = m->tam[day_idx];
        tair_pm = m->tpm[day_idx];
        sw_rad = m->par[day_idx] * PAR_2_SW;
        sw_rad_am = m->par_am[day_idx] * PAR_2_SW;
        sw_rad_pm = m->par_pm[day_idx] * PAR_2_SW;
        rain = m->rain[day_idx];
        vpd_am = m->vpd_am[day_idx] * KPA_2_PA;
        vpd_pm = m->vpd_pm[day_idx] * KPA_2_PA;
        wind_am = m->wind_am[day_idx];
        wind_pm = m->wind_pm[day_idx];
        wind_day = m->wind[day_idx];
        press = m->press[day_idx] * KPA_2_PA;
    }

    interception = calc_infiltration(p, s, rain);
    net_rad = calc_net_radiation(p, sw_rad, tair);
    soil_evap = calc_soil_evaporation(p, s, net_rad, press, tair);
    if (c->sub_daily) {
        soil_evap *= MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
    } else {
        net_rad_am = calc_net_radiation(p, sw_rad_am, tair);
        net_rad_pm = calc_net_radiation(p, sw_rad_pm, tair);
        soil_evap *= MOLE_WATER_2_G_WATER * G_TO_KG * (60.0 * 60.0 * daylen);
    }

    if (c->sub_daily) {
        /* mol m-2 s-1 to mm/30 min */
        transpiration = trans_leaf * MOLE_WATER_2_G_WATER * G_TO_KG * \
                        SEC_2_HLFHR;

    } else {
        /* umol m-2 s-1 */
        gpp_am = f->gpp_am * GRAMS_C_TO_MOL_C * MOL_TO_UMOL * HALF_DAY_2_SEC;
        gpp_pm = f->gpp_pm * GRAMS_C_TO_MOL_C * MOL_TO_UMOL * HALF_DAY_2_SEC;

        penman_canopy_wrapper(p, s, press, vpd_am, tair_am, wind_am, net_rad_am,
                              ca, gpp_am, &ga_am, &gs_am, &transpiration_am,
                              &LE_am, &omega_am);
        penman_canopy_wrapper(p, s, press, vpd_pm, tair_pm, wind_pm, net_rad_pm,
                              ca, gpp_pm, &ga_pm, &gs_pm, &transpiration_pm,
                              &LE_pm, &omega_pm);

        /* mol m-2 s-1 to mm/day */
        transpiration = (transpiration_am * MOLE_WATER_2_G_WATER * G_TO_KG * \
                         SEC_2_HALF_DAY) + \
                        (transpiration_pm * MOLE_WATER_2_G_WATER * G_TO_KG * \
                         SEC_2_HALF_DAY);

        f->omega = (omega_am + omega_pm) / 2.0;

        /* output in mol H20 m-2 s-1 */
        f->gs_mol_m2_sec = gs_am + gs_pm;
        f->ga_mol_m2_sec = ga_am + ga_pm;
    }

    et = transpiration + soil_evap + interception;
    update_water_storage(c, f, p, s, rain, interception, &transpiration,
                         &soil_evap, &et, &runoff);

    if (c->sub_daily) {
        sum_hourly_water_fluxes(f, soil_evap, transpiration, et,
                                interception, runoff, omega_leaf);
    } else {
        update_daily_water_struct(f, soil_evap, transpiration, et,
                                  interception, runoff);
    }

    return;
}

void update_water_storage(control *c, fluxes *f, params *p, state *s,
                          double rain, double interception,
                          double *transpiration, double *soil_evap,
                          double *et, double *runoff) {
    /* Calculate root and top soil plant available water and runoff.

    Soil drainage is estimated using a "leaky-bucket" approach with two
    soil layers. In reality this is a combined drainage and runoff
    calculation, i.e. "outflow". There is no drainage out of the "bucket"
    soil.

    Returns:
    --------
    outflow : float
        outflow [mm d-1]
    */
    double trans_frac, previous;

    /* reduce transpiration from the top soil if it is dry */
    trans_frac = p->fractup_soil * s->wtfac_topsoil;

    /* Total soil layer */
    s->pawater_topsoil += (rain - interception) - \
                          (*transpiration * trans_frac) - \
                           *soil_evap;

    if (s->pawater_topsoil < 0.0) {
        s->pawater_topsoil = 0.0;
    } else if (s->pawater_topsoil > p->wcapac_topsoil) {
        s->pawater_topsoil = p->wcapac_topsoil;
    }

    /* Total root zone */
    previous = s->pawater_root;
    s->pawater_root += (rain - interception) - *transpiration - *soil_evap;

    /* calculate runoff and remove any excess from rootzone */
    if (s->pawater_root > p->wcapac_root) {
        *runoff = s->pawater_root - p->wcapac_root;
        s->pawater_root -= *runoff;
    } else {
        *runoff = 0.0;
    }

    if (s->pawater_root < 0.0) {
        *transpiration = 0.0;
        *soil_evap = 0.0;
        *et = interception;
    }

    if (s->pawater_root < 0.0)
        s->pawater_root = 0.0;
    else if (s->pawater_root > p->wcapac_root)
        s->pawater_root = p->wcapac_root;

    s->delta_sw_store = s->pawater_root - previous;

    if (c->water_stress) {
        /* Calculate the soil moisture availability factors [0,1] in the
           topsoil and the entire root zone */
        calculate_soil_water_fac(c, p, s);
    } else {
        /* really this should only be a debugging option! */
        s->wtfac_topsoil = 1.0;
        s->wtfac_root = 1.0;
    }

    return;
}



double calc_infiltration(params *p, state* s, double rain) {
    /* Estimate "effective" rain, or infiltration I guess.

    Simple assumption that infiltration relates to leaf area
    and therefore canopy storage capacity (wetloss). Interception is
    likely to be ("more") erroneous if a canopy is subject to frequent daily
    rainfall I would suggest.

    Parameters:
    -------
    rain : float
        rainfall [mm d-1]

    */
    double interception;

    if (s->lai > 0.0) {
        interception = (rain * p->intercep_frac * \
                        MIN(1.0, s->lai / p->max_intercep_lai));
    } else {
        interception = 0.0;
    }

    return (interception);
}

double calc_soil_evaporation(params *p, state *s, double net_rad, double press,
                             double tair) {
    /* Use Penman eqn to calculate top soil evaporation flux at the
    potential rate.

    Soil evaporation is dependent upon soil wetness and plant cover. The net
    radiation term is scaled for the canopy cover passed to this func and
    the impact of soil wetness is accounted for in the wtfac term. As the
    soil dries the evaporation component reduces significantly.

    Key assumptions from Ritchie...

    * When plant provides shade for the soil surface, evaporation will not
    be the same as bare soil evaporation. Wind speed, net radiation and VPD
    will all belowered in proportion to the canopy density. Following
    Ritchie role ofwind, VPD are assumed to be negligible and are therefore
    ignored.

    These assumptions are based on work with crops and whether this holds
    for tree shading where the height from the soil to the base of the
    crown is larger is questionable.

    units = (mm/day)

    References:
    -----------
    * Ritchie, 1972, Water Resources Research, 8, 1204-1213.

    Parameters:
    -----------
    tair : float
        temperature [degC]
    net_rad : float
        net radiation [W m-2]
    press : float
        air pressure [kPa]

    Returns:
    --------
    soil_evap : float
        soil evaporation [mm d-1]

    */
    double lambda, gamma, slope, arg1, arg2, soil_evap;

    lambda = calc_latent_heat_of_vapourisation(tair);
    gamma = calc_pyschrometric_constant(press, lambda);
    slope = calc_slope_of_sat_vapour_pressure_curve(tair);

    /* mol H20 m-2 s-1 */
    soil_evap = ((slope / (slope + gamma)) * net_rad) / lambda;

    /*
      Surface radiation is reduced by overstory LAI cover. This empirical
      fit comes from Ritchie (1972) and is formed by a fit between the LAI
      of 5 crops types and the fraction of observed net radiation at the
      surface. Whilst the LAI does cover a large range, nominal 0–6, there
      are only 12 measurements and only three from LAI > 3. So this might
      not hold as well for a forest canopy?
      Ritchie 1972, Water Resources Research, 8, 1204-1213.
    */
    if (s->lai > 0.0)
        soil_evap *= exp(-0.398 * s->lai);

    /* reduce soil evaporation if top soil is dry */
    soil_evap *= s->wtfac_topsoil;

    return (soil_evap);
}

double calc_net_radiation(params *p, double sw_rad, double tair) {

    double net_rad, net_lw;

    /* Net loss of long-wave radn, Monteith & Unsworth '90, pg 52, eqn 4.17 */
    net_lw = 107.0 - 0.3 * tair;            /* W m-2 */

    /* Net radiation recieved by a surf, Monteith & Unsw '90, pg 54 eqn 4.21
        - note the minus net_lw is correct as eqn 4.17 is reversed in
          eqn 4.21, i.e Lu-Ld vs. Ld-Lu
        - NB: this formula only really holds for cloudless skies!
        - Bounding to zero, as we can't have negative soil evaporation, but you
          can have negative net radiation.
        - units: W m-2
    */
    net_rad = MAX(0.0, (1.0 - p->albedo) * sw_rad - net_lw);

    return (net_rad);
}


void penman_canopy_wrapper(params *p, state *s, double press, double vpd,
                           double tair, double wind, double rnet, double ca,
                           double gpp, double *ga, double *gsv,
                           double *transpiration, double *LE, double *omega) {
    /*
        Calculates transpiration at the canopy scale (or big leaf) using the
        Penman-Monteith

        Parameters:
        ----------
        parms : structure
            parameters
        state : structure
            state variables
        press : float
            atmospheric pressure (Pa)
        vpd : float
            vapour pressure deficit of air (Pa)
        tair : float
            air temperature (deg C)
        wind : float
            wind speed (m s-1)
        rnet : float
            net radiation (J m-2 s-1)
        ca : float
            ambient CO2 concentration (umol mol-1)
        gpp : float
            gross primary productivity (umol m-2 s-1)
        ga : float
            canopy scale boundary layer conductance (mol m-2 s-1; returned)
        gsv : float
            stomatal conductance to H20 (mol m-2 s-1; returned)
        transpiration : float
            transpiration (mol H2O m-2 s-1; returned)
        LE : float
            latent heat flux (W m-2; returned)


    */
    double gb, gv, gsc, epsilon, arg1, arg2, slope, gamma, lambda;

    /* stomtal conductance to CO2 */
    gsc = calc_stomatal_conductance(p, s, vpd, ca, gpp);

    /* stomtal conductance to H2O */
    *gsv = GSVGSC * gsc;

    *ga = canopy_boundary_layer_conduct(p, s->canht, wind, press, tair);
    gv = 1.0 / (1.0 / *gsv + 1.0 / *ga);

    lambda = calc_latent_heat_of_vapourisation(tair);
    gamma = calc_pyschrometric_constant(press, lambda);
    slope = calc_slope_of_sat_vapour_pressure_curve(tair);

    penman_monteith(press, vpd, rnet, slope, lambda, gamma, ga, &gv,
                    transpiration, LE);

    /* Calculate decoupling coefficient (McNaughton and Jarvis 1986) */
    epsilon = slope / gamma;
    *omega = (1.0 + epsilon) / (1.0 + epsilon + *ga / *gsv);

    return;
}

void penman_leaf_wrapper(params *p, state *s, double press, double vpd,
                         double tair, double tleaf, double wind, double rnet,
                         double gsc, double *transpiration, double *LE,
                         double *gbc, double *gh, double *gv, double *omega) {
    /*
        Calculates transpiration by leaves using the Penman-Monteith

        Parameters:
        ----------
        press : float
            atmospheric pressure (Pa)
        rnet : float
            net radiation (J m-2 s-1)
        vpd : float
            vapour pressure deficit of air (Pa)
        tair : float
            air temperature (deg C)
        transpiration : float
            transpiration (mol H2O m-2 s-1) (returned)
        LE : float
            latent heat flux, W m-2 (returned)


    */
    double slope, epsilon, lambda, arg1, arg2, gradn, gbhu, gbhf, gbh,
           gbv, gsv, gamma, Tdiff, sensible_heat, ema, Tk;

    /* Radiation conductance (mol m-2 s-1) */
    gradn = calc_radiation_conductance(tair);

    /* Boundary layer conductance for heat - single sided, forced
       convection (mol m-2 s-1) */
    gbhu = calc_bdn_layer_forced_conduct(tair, press, wind, p->leaf_width);

    /* Boundary layer conductance for heat - single sided, free convection */
    gbhf = calc_bdn_layer_free_conduct(tair, tleaf, press, p->leaf_width);

    /* Total boundary layer conductance for heat */
    gbh = gbhu + gbhf;

    /* Total conductance for heat - two-sided */
    *gh = 2.0 * (gbh + gradn);

    /* Total conductance for water vapour */
    gbv = GBVGBH * gbh;
    gsv = GSVGSC * gsc;
    *gv = (gbv * gsv) / (gbv + gsv);
    *gbc = gbh / GBHGBC;

    lambda = calc_latent_heat_of_vapourisation(tair);
    gamma = calc_pyschrometric_constant(press, lambda);
    slope = calc_slope_of_sat_vapour_pressure_curve(tair);

    penman_monteith(press, vpd, rnet, slope, lambda, gamma, gh, gv,
                    transpiration, LE);

    /* Calculate decoupling coefficient (McNaughton and Jarvis 1986) */
    epsilon = slope / gamma;
    *omega = (1.0 + epsilon) / (1.0 + epsilon + gbv / gsv);

    return;
}

void penman_monteith(double press, double vpd, double rnet, double slope,
                     double lambda, double gamma, double *gh, double *gv,
                     double *transpiration, double *LE) {
    /*
        Calculates transpiration using the Penman-Monteith

        Parameters:
        ----------
        press : float
            atmospheric pressure (Pa)
        vpd : float
            vapour pressure deficit of air (Pa)
        rnet : float
            net radiation (J m-2 s-1)
        slope : float
            slope of VPD/T curve, Pa K-1
        lambda : flot
            latent heat of water at air T, J mol-1
        gamma : float
            psychrometric constant, J mol-1
        gh : float
            boundary layer conductance to heat (free & forced & radiative
            components), mol m-2 s-1
        gv : float
            conductance to water vapour (stomatal & bdry layer components),
            mol m-2 s-1
        transpiration : float
            transpiration (mol H2O m-2 s-1; returned)
        LE : float
            latent heat flux (W m-2; returned)
        omega : float
            decoupling coefficient (unitless; returned)
    */
    double arg1, arg2, epsilon;

    if (*gv > 0.0) {
        arg1 = slope * rnet + vpd * *gh * CP * MASS_AIR;
        arg2 = slope + gamma * *gh / *gv;
        *LE = arg1 / arg2; /* W m-2 */
        *transpiration = *LE / lambda; /* mol H20 m-2 s-1 */
    } else {
        *transpiration = 0.0;
    }

    /* Should not be negative - not sure gv>0.0 catches it as g0 = 1E-09? */
    *transpiration = MAX(0.0, *transpiration);

    return;
}

double calc_stomatal_conductance(params *p, state *s, double vpd, double Ca,
                                 double gpp) {
    /*
        Calculate stomatal conductance using Belinda's model. For the medlyn
        model this is already in conductance to CO2, so the 1.6 from the
        corrigendum to Medlyn et al 2011 is missing here

    References:
    -----------
    For conversion factor for conductance see...
    * Jones (1992) Plants and microclimate, pg 56 + Appendix 3
    * Diaz et al (2007) Forest Ecology and Management, 244, 32-40.

    Stomatal Model:
    * Medlyn et al. (2011) Global Change Biology, 17, 2134-2144.
    **Note** Corrigendum -> Global Change Biology, 18, 3476.

    Parameters:
    -----------
    g1 : float
        slope
    wtfac : float
        water availability scaler [0,1]
    vpd : float
        vapour pressure deficit (Pa)
    Ca : float
        atmospheric co2 [umol mol-1]
    gpp : float
        photosynthesis at the canopy scale (umol m-2 s-1)

    Returns:
    --------
    gs : float
        stomatal conductance (mol CO2 m-2 s-1)
    */
    double gs_over_a, g1, g0, gsc;

    g1 = p->g1 * s->wtfac_root;
    g0 = 0.0; /* p->g0; */
    gs_over_a = (1.0 + g1 / sqrt(vpd)) / Ca;
    gsc = MAX(g0, g0 + gs_over_a * gpp);

    /* mol m-2 s-1 */
    return (gsc);

}

double canopy_boundary_layer_conductance(params *p, double wind, double canht) {
    /*  Canopy boundary layer conductance, ga or 1/ra

    Characterises the heat/water vapour from evaporating surface, but does
    not account for leaf boundary layer conductance, which is the parellel
    sum of single leaf boundary layer conductances for all leaves in the
    canopy.

    Notes:
    ------
    'Estimates of ga for pine canopies from LAI of 3 to 6 vary from
    3.5 to 1.1 mol m-2 s-1  (Kelliher et al., 1993; Juang et al., 2007).'
    Drake et al, 2010, 17, pg. 1526.

    References:
    ------------
    * Jones 1992, pg. 67-8.
    * Monteith and Unsworth (1990), pg. 248. Note this in the inverted form
      of what is in Monteith (ga = 1 / ra)
    * Allen et al. (1989) pg. 651.
    * Gash et al. (1999) Ag forest met, 94, 149-158.

    Parameters:
    -----------
    wind : float
        average daytime wind speed [m s-1]

    Returns:
    --------
    ga : float
        canopy boundary layer conductance [m s-1]
    */

    /* z0m roughness length governing momentum transfer [m] */
    double z0m, z0h, d, arg1, arg2, arg3;
    double vk = 0.41;
    z0m = p->dz0v_dh * canht;

    /*
       z0h roughness length governing transfer of heat and vapour [m]
      *Heat tranfer typically less efficent than momentum transfer. There is
       a lot of variability in values quoted for the ratio of these two...
       JULES uses 0.1, Campbell and Norman '98 say z0h = z0m / 5. Garratt
       and Hicks, 1973/ Stewart et al '94 say z0h = z0m / 7. Therefore for
       the default I am following Monteith and Unsworth, by setting the
       ratio to be 1, the code below is identical to that on page 249,
       eqn 15.7
    */
    z0h = p->z0h_z0m * z0m;

    /* zero plan displacement height [m] */
    d = p->displace_ratio * canht;

    arg1 = (vk * vk) * wind;
    arg2 = log((canht - d) / z0m);
    arg3 = log((canht - d) / z0h);

    return (arg1 / (arg2 * arg3));
}

double calc_radiation_conductance(double tair) {
    /*  Returns the 'radiation conductance' at given temperature.

        Units: mol m-2 s-1

        References:
        -----------
        * Formula from Ying-Ping's version of Maestro, cf. Wang and Leuning
          1998, Table 1,
        * See also Jones (1992) p. 108.
        * And documented in Medlyn 2007, equation A3, although I think there
          is a mistake. It should be Tk**3 not Tk**4, see W & L.
    */
    double grad;
    double Tk;

    Tk = tair + DEG_TO_KELVIN;
    grad = 4.0 * SIGMA * (Tk * Tk * Tk) * LEAF_EMISSIVITY / (CP * MASS_AIR);

    return (grad);
}

double calc_bdn_layer_forced_conduct(double tair, double press, double wind,
                                     double leaf_width) {
    /*
        Boundary layer conductance for heat - single sided, forced convection
        (mol m-2 s-1)
        See Leuning et al (1995) PC&E 18:1183-1200 Eqn E1
    */
    double cmolar, Tk, gbh;

    Tk = tair + DEG_TO_KELVIN;
    cmolar = press / (RGAS * Tk);
    gbh = 0.003 * sqrt(wind / leaf_width) * cmolar;

    return (gbh);
}

double calc_bdn_layer_free_conduct(double tair, double tleaf, double press,
                                   double leaf_width) {
    /*
        Boundary layer conductance for heat - single sided, free convection
        (mol m-2 s-1)
        See Leuning et al (1995) PC&E 18:1183-1200 Eqns E3 & E4
    */
    double cmolar, Tk, gbh, grashof, leaf_width_cubed;
    double tolerance = 1E-08;

    Tk = tair + DEG_TO_KELVIN;
    cmolar = press / (RGAS * Tk);
    leaf_width_cubed = leaf_width * leaf_width * leaf_width;

    if (float_eq((tleaf - tair), 0.0)) {
        gbh = 0.0;
    } else {
        grashof = 1.6E8 * fabs(tleaf - tair) * leaf_width_cubed;
        gbh = 0.5 * DHEAT * pow(grashof, 0.25) / leaf_width * cmolar;
    }

    return (gbh);
}


double canopy_boundary_layer_conduct(params *p, double canht, double wind,
                                     double press, double tair) {
    /*  Canopy boundary layer conductance, ga (from Jones 1992 p 68)


    Notes:
    ------
    'Estimates of ga for pine canopies from LAI of 3 to 6 vary from
    3.5 to 1.1 mol m-2 s-1  (Kelliher et al., 1993; Juang et al., 2007).'
    Drake et al, 2010, 17, pg. 1526.

    References:
    ------------
    * Jones 1992, pg. 67-8.
    * Monteith and Unsworth (1990), pg. 248. Note this in the inverted form
      of what is in Monteith (ga = 1 / ra)
    * Allen et al. (1989) pg. 651.
    * Gash et al. (1999) Ag forest met, 94, 149-158.

    Parameters:
    -----------
    params : p
        parameters structure
    canht : float
        canopy height (m)
    wind : float
        wind speed (m s-1)
    press : float
        atmospheric pressure (Pa)
    tair : float
        air temperature (deg C)

    Returns:
    --------
    ga : float
        canopy boundary layer conductance (mol m-2 s-1)
    */

    /* z0m roughness length governing momentum transfer [m] */
    double z0m, z0h, d, arg1, arg2, arg3, tk, ga, cmolar;
    double vk = 0.41;

    tk = tair + DEG_TO_KELVIN;

    /* Convert from mm s-1 to mol m-2 s-1 */
    cmolar = press / (RGAS * tk);

    /* roughness length for momentum */
    z0m = p->dz0v_dh * canht;

    /*
       z0h roughness length governing transfer of heat and vapour [m]
      *Heat tranfer typically less efficent than momentum transfer. There is
       a lot of variability in values quoted for the ratio of these two...
       JULES uses 0.1, Campbell and Norman '98 say z0h = z0m / 5. Garratt
       and Hicks, 1973/ Stewart et al '94 say z0h = z0m / 7. Therefore for
       the default I am following Monteith and Unsworth, by setting the
       ratio to be 1, the code below is identical to that on page 249,
       eqn 15.7
    */
    z0h = p->z0h_z0m * z0m;

    /* zero plan displacement height [m] */
    d = p->displace_ratio * canht;

    arg1 = (vk * vk) * wind;
    arg2 = log((canht - d) / z0m);
    arg3 = log((canht - d) / z0h);

    ga = (arg1 / (arg2 * arg3)) * cmolar;

    return (ga);
}

double calc_slope_of_sat_vapour_pressure_curve(double tair) {
    /*
        Constant slope in Penman-Monteith equation

        Parameters:
        -----------
        tavg : float
            average daytime temperature

        Returns:
        --------
        slope : float
            slope of saturation vapour pressure curve [Pa K-1]

    */
    double arg1, arg2, slope;

    /* Const slope in Penman-Monteith equation  (Pa K-1) */
    arg1 = calc_sat_water_vapour_press(tair + 0.1);
    arg2 = calc_sat_water_vapour_press(tair);
    slope = (arg1 - arg2) / 0.1;

    return (slope);
}


double calc_pyschrometric_constant(double press, double lambda) {
    /* Psychrometric constant ratio of specific heat of moist air at
    a constant pressure to latent heat of vaporisation.

    Parameters:
    -----------
    press : float
        air pressure (Pa)
    lambda : float
         latent heat of water vaporization (J mol-1)


    Returns:
    --------
    gamma : float
        pyschrometric_constant [J mol-1]

    */
    return ( CP * MASS_AIR * press / lambda );

}

double calc_latent_heat_of_vapourisation(double tair) {
    /*   Latent heat of water vapour at air temperature

        Returns:
        -----------
        lambda : float
            latent heat of water vaporization [J mol-1]
    */
    return ( (H2OLV0 - 2.365E3 * tair) * H2OMW );

}

double calc_sat_water_vapour_press(double tac) {
    /*
        Calculate saturated water vapour pressure (Pa) at
        temperature TAC (Celsius). From Jones 1992 p 110 (note error in
        a - wrong units)
    */
    return (613.75 * exp(17.502 * tac / (240.97 + tac)));
}

void initialise_soil_moisture_parameters(control *c, params *p) {
    /*
      initialise parameters, if these are not known for the site use
      values derived from Cosby et al to calculate the amount of plant
      available water.
     */

    double theta_fc_topsoil, theta_wp_topsoil, theta_fc_root, theta_wp_root;
    double *fsoil_top = NULL, *fsoil_root = NULL;

    if (c->calc_sw_params) {
        fsoil_top = get_soil_fracs(p->topsoil_type);
        fsoil_root = get_soil_fracs(p->rootsoil_type);

        /* top soil */
        calc_soil_params(fsoil_top, &theta_fc_topsoil, &theta_wp_topsoil,
                         &p->theta_sat_topsoil, &p->b_topsoil, &p->psi_sat_topsoil);

        /* Plant available water in top soil (mm) */
        p->wcapac_topsoil = p->topsoil_depth * (theta_fc_topsoil - theta_wp_topsoil);

        /* Root zone */
        calc_soil_params(fsoil_root, &theta_fc_root, &theta_wp_root,
                         &p->theta_sat_root, &p->b_root, &p->psi_sat_root);

        /* Plant available water in rooting zone (mm) */
        p->wcapac_root = p->rooting_depth * (theta_fc_root - theta_wp_root);
    }

    /* calculate Landsberg and Waring SW modifier parameters if not
       specified by the user based on a site calibration */
    if (p->ctheta_topsoil < -900.0 && p->ntheta_topsoil  < -900.0 &&
        p->ctheta_root < -900.0 && p->ntheta_root < -900.0) {
        get_soil_params(p->topsoil_type, &p->ctheta_topsoil, &p->ntheta_topsoil);
        get_soil_params(p->rootsoil_type, &p->ctheta_root, &p->ntheta_root);
    }
    /*
    printf("%f\n", p->wcapac_topsoil);
    printf("%f\n\n", p->wcapac_root);

    printf("%f\n", p->ctheta_topsoil);
    printf("%f\n", p->ntheta_topsoil);
    printf("%f\n", p->ctheta_root);
    printf("%f\n", p->ntheta_root);
    printf("%f\n", p->rooting_depth);

    exit(1); */



    free(fsoil_top);
    free(fsoil_root);

    return;
}


double *get_soil_fracs(char *soil_type) {
    /* Based on Table 2 in Cosby et al 1984, page 2.
    Fractions of silt, sand and clay (in that order)
    */
    double *fsoil = malloc(3 * sizeof(double));

    if (strcmp(soil_type, "sand") == 0) {
        fsoil[0] = 0.05;
        fsoil[1] = 0.92;
        fsoil[2] = 0.03;
    } else if (strcmp(soil_type, "loamy_sand") == 0) {
        fsoil[0] = 0.12;
        fsoil[1] = 0.82;
        fsoil[2] = 0.06;
    } else if (strcmp(soil_type, "sandy_loam") == 0) {
        fsoil[0] = 0.32;
        fsoil[1] = 0.58;
        fsoil[2] = 0.1;
    } else if (strcmp(soil_type, "loam") == 0) {
        fsoil[0] = 0.39;
        fsoil[1] = 0.43;
        fsoil[2] = 0.18;
    } else if (strcmp(soil_type, "silty_loam") == 0) {
        fsoil[0] = 0.7;
        fsoil[1] = 0.17;
        fsoil[2] = 0.13;
    } else if (strcmp(soil_type, "sandy_clay_loam") == 0) {
        fsoil[0] = 0.15;
        fsoil[1] = 0.58;
        fsoil[2] = 0.27;
    } else if (strcmp(soil_type, "clay_loam") == 0) {
        fsoil[0] = 0.34;
        fsoil[1] = 0.32;
        fsoil[2] = 0.34;
    } else if (strcmp(soil_type, "silty_clay_loam") == 0) {
        fsoil[0] = 0.56;
        fsoil[1] = 0.1;
        fsoil[2] = 0.34;
    } else if (strcmp(soil_type, "sandy_clay") == 0) {
        fsoil[0] = 0.06;
        fsoil[1] = 0.52;
        fsoil[2] = 0.42;
    } else if (strcmp(soil_type, "silty_clay") == 0) {
        fsoil[0] = 0.47;
        fsoil[1] = 0.06;
        fsoil[2] = 0.47;
    } else if (strcmp(soil_type, "clay") == 0) {
        fsoil[0] = 0.2;
        fsoil[1] = 0.22;
        fsoil[2] = 0.58;
    } else {
        prog_error("Could not understand soil type", __LINE__);
    }

    return (fsoil);
}

void get_soil_params(char *soil_type, double *c_theta, double *n_theta) {
    /* For a given soil type, get the parameters for the soil
    moisture availability based on Landsberg and Waring, with updated
    parameters from Landsberg and Sands (2011), pg 190, Table 7.1

    Table also has values from Saxton for soil texture, perhaps makes more
    sense to use those than Cosby? Investigate?

    Reference
    ---------
    * Landsberg and Sands (2011) Physiological ecology of forest production.
    * Landsberg and Waring (1997) Forest Ecology & Management, 95, 209-228.
    */

    if (strcmp(soil_type, "clay") == 0) {
        *c_theta = 0.4;
        *n_theta = 3.0;
    } else if (strcmp(soil_type, "clay_loam") == 0) {
        *c_theta = 0.5;
        *n_theta = 5.0;
    } else if (strcmp(soil_type, "loam") == 0) {
        *c_theta = 0.55;
        *n_theta = 6.0;
    } else if (strcmp(soil_type, "loamy_sand") == 0) {
        *c_theta = 0.65;
        *n_theta = 8.0;
    } else if (strcmp(soil_type, "sand") == 0) {
        *c_theta = 0.7;
        *n_theta = 9.0;
    } else if (strcmp(soil_type, "sandy_clay") == 0) {
        *c_theta = 0.45;
        *n_theta = 4.0;
    } else if (strcmp(soil_type, "sandy_clay_loam") == 0) {
        *c_theta = 0.525;
        *n_theta = 5.5;
    } else if (strcmp(soil_type, "sandy_loam") == 0) {
        *c_theta = 0.6;
        *n_theta = 7.0;
    } else if (strcmp(soil_type, "silt") == 0) {
        *c_theta = 0.625;
        *n_theta = 7.5;
    } else if (strcmp(soil_type, "silty_clay") == 0) {
        *c_theta = 0.425;
        *n_theta = 3.5;
    } else if (strcmp(soil_type, "silty_clay_loam") == 0) {
        *c_theta = 0.475;
        *n_theta = 4.5;
    } else if (strcmp(soil_type, "silty_loam") == 0) {
        *c_theta = 0.575;
        *n_theta = 6.5;
    } else {
        prog_error("There are no parameters for your soil type", __LINE__);
    }

    return;
}

void calc_soil_params(double *fsoil, double *theta_fc, double *theta_wp,
                      double *theta_sp, double *b, double *psi_sat_mpa) {
    /* Cosby parameters for use within the Clapp Hornberger soil hydraulics
    scheme are calculated based on the texture components of the soil.

    NB: Cosby et al were ambiguous in their paper as to what log base to
    use.  The correct implementation is base 10, as below.

    Parameters:
    ----------
    fsoil : list
        fraction of silt, sand, and clay (in that order

    Returns:
    --------
    theta_fc : float
        volumetric soil water concentration at field capacity
    theta_wp : float
        volumetric soil water concentration at the wilting point

    */
    /* soil suction of 3.364m and 152.9m, or equivalent of -0.033 & -1.5 MPa */
    double pressure_head_wilt = -152.9;
    double pressure_head_crit = -3.364;
    double KPA_2_MPA, METER_OF_HEAD_TO_MPA, psi_sat;

    /* *Note* subtle unit change to be consistent with fractions as opposed
      to percentages of sand, silt, clay, e.g. I've changed the slope in
      the "b" Clapp paramter from 0.157 to 15.7

      Also Cosby is unclear about which log base were used. 'Generally' now
      assumed that logarithms to the base 10
    */

    /* Clapp Hornberger exponent [-] */
    *b = 3.1 + 15.7 * fsoil[CLAY] - 0.3 * fsoil[SAND];

    /* soil matric potential at saturation, taking inverse of log (base10)
      units = m (0.01 converts from mm to m) */
    psi_sat = 0.01 * -(pow(10.0, (1.54 - 0.95 * fsoil[SAND] + 0.63 * fsoil[SILT])));

    /* Height (m) x gravity (m/s2) = pressure (kPa) */
    KPA_2_MPA = 0.001;
    METER_OF_HEAD_TO_MPA = 9.81 * KPA_2_MPA;
    *psi_sat_mpa = psi_sat * METER_OF_HEAD_TO_MPA;

    /* volumetric soil moisture concentrations at the saturation point */
    *theta_sp = 0.505 - 0.037 * fsoil[CLAY] - 0.142 * fsoil[SAND];

    /* volumetric soil moisture concentrations at the wilting point
       assumed to equal suction of -1.5 MPa or a depth of water of 152.9 m */
    *theta_wp = *theta_sp * pow((psi_sat / pressure_head_wilt), (1.0 / *b));

    /* volumetric soil moisture concentrations at field capacity assumed to
       equal a suction of -0.0033 MPa or a depth of water of 3.364 m */
    *theta_fc = *theta_sp * pow((psi_sat / pressure_head_crit), (1.0 / *b));

    return;

}

void calculate_soil_water_fac(control *c, params *p, state *s) {
    /* Estimate a relative water availability factor [0..1]

    A drying soil results in physiological stress that can induce stomatal
    closure and reduce transpiration. Further, N mineralisation depends on
    top soil moisture.

    s->qs = 0.2 in SDGVM

    References:
    -----------
    * Landsberg and Waring (1997) Forest Ecology and Management, 95, 209-228.
      See --> Figure 2.
    * Egea et al. (2011) Agricultural Forest Meteorology, 151, 1370-1384.

    But similarly see:
    * van Genuchten (1981) Soil Sci. Soc. Am. J, 44, 892--898.
    * Wang and Leuning (1998) Ag Forest Met, 91, 89-111.

    * Pepper et al. (2008) Functional Change Biology, 35, 493-508

    Returns:
    --------
    wtfac_topsoil : float
        water availability factor for the top soil [0,1]
    wtfac_root : float
        water availability factor for the root zone [0,1]
    */

    double smc_topsoil, smc_root, psi_swp_topsoil, arg1, arg2, arg3,
           psi_swp_root, b;

    /* turn into fraction... */
    smc_topsoil = s->pawater_topsoil / p->wcapac_topsoil;
    smc_root = s->pawater_root / p->wcapac_root;

    if (c->sw_stress_model == 0) {
        s->wtfac_topsoil = pow(smc_topsoil, p->qs);
        s->wtfac_root = pow(smc_root, p->qs);

    } else if (c->sw_stress_model == 1) {
        s->wtfac_topsoil = calc_sw_modifier(smc_topsoil, p->ctheta_topsoil,
                                            p->ntheta_topsoil);


        s->wtfac_root = calc_sw_modifier(smc_root, p->ctheta_root,
                                         p->ntheta_root);

    } else if (c->sw_stress_model == 2) {

        /* Stomatal limitaiton
           Exponetial function to reduce g1 with soil water limitation
           based on Zhou et al. 2013, AFM, following Makela et al 1996.
           For the moment I have hardwired the PFT parameter as I am still
           testing.
           Because the model is a daily model we are assuming that LWP is
           well approximated by the night SWP. */

        if (float_eq(smc_topsoil, 0.0)) {
            psi_swp_topsoil = -1.5;
        } else {
            arg1 = s->psi_sat_topsoil;
            arg2 = smc_topsoil /s->theta_sat_topsoil;
            arg3 = s->b_topsoil;
            psi_swp_topsoil = arg1 * pow(arg2, arg3);
        }

        if (float_eq(smc_root, 0.0)) {
            psi_swp_root = -1.5;
        } else {
            arg1 = s->psi_sat_root;
            arg2 = smc_root/s->theta_sat_root;
            arg3 = s->b_root;
            psi_swp_root = arg1 * pow(arg2, arg3);
        }

        /* multipliy these by g1, same as eqn 3 in Zhou et al. 2013. */
        b = 0.66;

        s->wtfac_topsoil = exp(b * psi_swp_topsoil);
        s->wtfac_root = exp(b * psi_swp_root);
    }

    return;
}

double calc_sw_modifier(double theta, double c_theta, double n_theta) {
    /* From Landsberg and Waring */
    return (1.0  / (1.0 + pow(((1.0 - theta) / c_theta), n_theta)));
}

void sum_hourly_water_fluxes(fluxes *f, double soil_evap_hlf_hr,
                             double transpiration_hlf_hr, double et_hlf_hr,
                             double interception_hlf_hr,
                             double runoff_hlf_hr, double omega_hlf_hr) {

    /* add half hour fluxes to day total store */
    f->soil_evap += soil_evap_hlf_hr;
    f->transpiration += transpiration_hlf_hr;
    f->et += et_hlf_hr;
    f->interception += interception_hlf_hr;
    f->runoff += runoff_hlf_hr;
    f->omega += omega_hlf_hr; /* average at the end of hour loop */

    return;
}

void update_daily_water_struct(fluxes *f, double day_soil_evap,
                               double day_transpiration, double day_et,
                               double day_interception, double day_runoff) {

    /* add half hour fluxes to day total store */
    f->soil_evap = day_soil_evap;
    f->transpiration = day_transpiration;
    f->et = day_et;
    f->interception = day_interception;
    f->runoff = day_runoff;

    return;
}

void zero_water_day_fluxes(fluxes *f) {

    f->et = 0.0;
    f->soil_evap = 0.0;
    f->transpiration = 0.0;
    f->interception = 0.0;
    f->runoff = 0.0;
    f->gs_mol_m2_sec = 0.0;

    return;
}
