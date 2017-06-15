#include "water_balance_sub_daily.h"

void initialise_soils_sub_daily(control *c, fluxes *f, params *p, state *s) {
    /*
        Initialise soil water state & parameters
            - We have two options: (i) simple two-layer bucket, or
            - multi-layer soil model.

        Currently, I'm just using the same soil type for all layers in the
        multi-layer model, once it all works we can add something to allow
        texture to be changes by layer/depth.

        NB. get_soil_fracs, calc_soil_params, get_soil_params live in
            water_balance.c
    */

    double *fsoil_top = NULL, *fsoil_root = NULL;
    int     i;

    /* site params not known, so derive them based on Cosby et al */
    if (c->calc_sw_params) {
        fsoil_top = get_soil_fracs(p->topsoil_type);
        fsoil_root = get_soil_fracs(p->rootsoil_type);

        /* top soil */
        calc_soil_params(fsoil_top, &p->theta_fc_topsoil, &p->theta_wp_topsoil,
                         &p->theta_sp_topsoil, &p->b_topsoil,
                         &p->psi_sat_topsoil);

        /* Plant available water in top soil (mm) */
        p->wcapac_topsoil = p->topsoil_depth  *\
                            (p->theta_fc_topsoil - p->theta_wp_topsoil);
        /* Root zone */
        calc_soil_params(fsoil_root, &p->theta_fc_root, &p->theta_wp_root,
                         &p->theta_sp_root, &p->b_root, &p->psi_sat_root);

        /* Plant available water in rooting zone (mm) */
        p->wcapac_root = p->rooting_depth * \
                            (p->theta_fc_root - p->theta_wp_root);
    }


    /* calculate Landsberg and Waring SW modifier parameters if not
       specified by the user based on a site calibration */
    if (p->ctheta_topsoil < -900.0 && p->ntheta_topsoil  < -900.0 &&
        p->ctheta_root < -900.0 && p->ntheta_root < -900.0) {
        get_soil_params(p->topsoil_type, &p->ctheta_topsoil, &p->ntheta_topsoil);
        get_soil_params(p->rootsoil_type, &p->ctheta_root, &p->ntheta_root);
    }

    /* Set up all the hydraulics stuff */
    if (c->water_balance == HYDRAULICS) {
        calc_saxton_stuff(p, fsoil_root);

        for (i = 0; i < p->wetting; i++) {
            s->wetting_bot[i] = 0.0;
            s->wetting_top[i] = 0.0;
        }

        /* saturate the top layer */
        s->wetting_bot[0] = s->thickness[0];

        /* Initalise SW fraction - we should read this from param file */
        s->initial_water = 0.0;
        for (i = 0; i < p->soil_layers; i++) {
            s->water_frac[i] = 0.4;
            s->initial_water += 1E3 * (s->water_frac[i] * s->thickness[i]);
        }

        /*
        ** The loop needs to be outside the func as we need to be able to
        ** calculate the soil conductance per layer and call this via
        ** the integration func when we update the soil water balance
        */
        for (i = 0; i < p->core; i++) {
            f->soil_conduct[i] = calc_soil_conductivity(s->water_frac[i],
                                                        p->cond1[i], p->cond2[i],
                                                        p->cond3[i]);
        }

        calc_soil_root_resistance(f, p, s);
        calc_soil_water_potential(f, p, s);

        /* Calculate the weighted soil-water-potential */
        calc_water_uptake_per_layer(f, p, s);
    }

    if (c->calc_sw_params) {
        free(fsoil_top);
        free(fsoil_root);
    }

    return;
}

void calculate_water_balance_sub_daily(control *c, canopy_wk *cw, fluxes *f,
                                       met *m, nrutil *nr, params *p, state *s,
                                       int daylen, double trans,
                                       double omega_leaf, double rnet_leaf,
                                       double et_deficit, double year,
                                       double doy) {
    /*
        Calculate the water balance (including all water fluxes).
        - we are using all the hydraulics instead

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
            length of day in hours. (Dummy argument, only passed for daily
            model)
        trans : double
            transpiration (Dummy argument, only passed for sub-daily
            model)
        omega_leaf : double
            decoupling coefficient (Dummy argument, only passed for sub-daily
            model)
        rnet_leaf : double
            total canopy rnet (Dummy argument, only passed for sub-daily model)
    */

    int    i;
    double soil_evap, et, interception, runoff, conv, transpiration, net_rad;
    double SEC_2_DAY, DAY_2_SEC, transpiration_am, transpiration_pm, gs_am;
    double canopy_evap, surface_water;

    // Water drained through the bottom soil layer
    double water_lost = 0.0;

    if (c->water_balance == HYDRAULICS) {

        zero_water_movement(f, p);

        // calculate potential canopy evap rate, this may be reduced later
        // depending on canopy water storage
        canopy_evap = calc_canopy_evaporation(m, p, s, rnet_leaf);

        /* mol m-2 s-1 to mm d-1 */
        conv = MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
        canopy_evap *= conv;

        /* We could now replace this interception bit with the Rutter scheme? */
        calc_interception(c, m, p, f, s, &surface_water, &interception,
                          &canopy_evap);

        net_rad = calc_net_radiation(p, m->sw_rad, m->tair);
        soil_evap = calc_soil_evaporation(m, p, s, net_rad);
        soil_evap *= MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;

        /* mol m-2 s-1 to mm/30 min */
        transpiration = trans * MOLE_WATER_2_G_WATER * G_TO_KG * \
                        SEC_2_HLFHR;

        et = transpiration + soil_evap + canopy_evap;

        //
        // The loop needs to be outside the func as we need to be able to
        // calculate the soil conductance per layer and call this via
        // the integration func when we update the soil water balance
        //
        for (i = 0; i < p->core; i++) {
            f->soil_conduct[i] = calc_soil_conductivity(s->water_frac[i],
                                                        p->cond1[i],
                                                        p->cond2[i],
                                                        p->cond3[i]);
        }

        calc_soil_water_potential(f, p, s);
        calc_soil_root_resistance(f, p, s);

        // If we have leaves we are transpiring
        if (s->lai > 0.0) {
            calc_water_uptake_per_layer(f, p, s);
        }

        // Calculates the thickness of the top dry layer and determines water
        // lost in upper layers due to evaporation
        calc_wetting_layers(f, p, s, soil_evap, surface_water);
        extract_water_from_layers(f, s, soil_evap, transpiration);

        //
        // determines water movement between soil layers due drainage
        // down the profile
        //
        for (i = 0; i < p->soil_layers; i++) {
            if (c->soil_drainage == GRAVITY) {
                calc_soil_balance(f, nr, p, s, i, &water_lost);
            } else if (c->soil_drainage == CASCADING) {
                // Redistribute soil water following a cascading or
                // 'tipping bucket' approach, much simpler and computational
                // effective. We have made an assumption about the drainage
                // rate to make this work
                calc_soil_balance_cascading(f, nr, p, s, i, &water_lost);
            }
        }

        //
        // how much surface water infiltrantes the first soil layer in the
        // current time step? Water which does not infiltrate in a single step
        // is considered runoff
        //
        runoff = calc_infiltration(f, p, s, surface_water);
        runoff += water_lost * M_TO_MM;

        // Don't see point of calculating these again
        // Find SWP & soil resistance without updating waterfrac yet
        calc_soil_water_potential(f, p, s);
        calc_soil_root_resistance(f, p, s);

        update_soil_water_storage(f, p, s, &soil_evap, &transpiration);
        et = transpiration + soil_evap + canopy_evap;
    } else {

        // Simple soil water bucket appoximation

        // calculate potential canopy evap rate, this may be reduced later
        // depending on canopy water storage
        canopy_evap = calc_canopy_evaporation(m, p, s, rnet_leaf);

        /* mol m-2 s-1 to mm/day */
        conv = MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
        canopy_evap *= conv;
        calc_interception(c, m, p, f, s, &surface_water, &interception,
                          &canopy_evap);

        net_rad = calc_net_radiation(p, m->sw_rad, m->tair);
        soil_evap = calc_soil_evaporation(m, p, s, net_rad);
        soil_evap *= MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;

        /* mol m-2 s-1 to mm/30 min */
        transpiration = trans * MOLE_WATER_2_G_WATER * G_TO_KG * \
                        SEC_2_HLFHR;

        /*
        ** NB. et, transpiration & soil evap may all be adjusted in
        ** update_water_storage if we don't have sufficient water
        */
        et = transpiration + soil_evap + canopy_evap;


        update_water_storage(c, f, p, s, surface_water, interception, canopy_evap,
                             &transpiration, &soil_evap, &et, &runoff);

    }

    if (c->water_store) {
        // Do we need to take any water from the plant store? This function
        // also checks to for drought-induced mortality
        update_plant_water_store(cw, p, s, &transpiration, &et, et_deficit,
                                 year, doy);
    }

    sum_hourly_water_fluxes(f, soil_evap, transpiration, et, interception,
                            surface_water, canopy_evap, runoff, omega_leaf,
                            m->rain);

}

void zero_water_movement(fluxes *f, params *p) {

    // Losses and gains of water both from PPT and between layers need to be
    // zero'd at the start of each timestep
    int i;

    for (i = 0; i < p->core; i++) {
        f->water_loss[i] = 0.0;
        f->water_gain[i] = 0.0;
        f->ppt_gain[i] = 0.0;
        f->fraction_uptake[i] = 0.0;
        f->est_evap[i] = 0.0;
    }

    return;
}

void setup_hydraulics_arrays(fluxes *f, params *p, state *s) {
    /* Allocate the necessary memory for all the hydraulics arrays */
    p->potA = malloc(p->core * sizeof(double));
    if (p->potA == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's potA\n");
        exit(EXIT_FAILURE);
    }

    p->potB = malloc(p->core * sizeof(double));
    if (p->potB == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's potB\n");
        exit(EXIT_FAILURE);
    }

    p->cond1 = malloc(p->core * sizeof(double));
    if (p->cond1 == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's cond1\n");
        exit(EXIT_FAILURE);
    }

    p->cond2 = malloc(p->core * sizeof(double));
    if (p->cond1 == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's cond2\n");
        exit(EXIT_FAILURE);
    }

    p->cond3 = malloc(p->core * sizeof(double));
    if (p->cond1 == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's cond3\n");
        exit(EXIT_FAILURE);
    }

    p->porosity = malloc(p->core * sizeof(double));
    if (p->porosity == NULL) {
        fprintf(stderr, "malloc failed allocating porosity\n");
        exit(EXIT_FAILURE);
    }

    p->field_capacity = malloc(p->core * sizeof(double));
    if (p->field_capacity == NULL) {
        fprintf(stderr, "malloc failed allocating field_capacity\n");
        exit(EXIT_FAILURE);
    }

    f->soil_conduct = malloc(p->core * sizeof(double));
    if (f->soil_conduct == NULL) {
        fprintf(stderr, "malloc failed allocating soil_conduct\n");
        exit(EXIT_FAILURE);
    }

    f->swp = malloc(p->core * sizeof(double));
    if (f->swp == NULL) {
        fprintf(stderr, "malloc failed allocating swp\n");
        exit(EXIT_FAILURE);
    }

    f->soilR = malloc(p->core * sizeof(double));
    if (f->soilR == NULL) {
        fprintf(stderr, "malloc failed allocating soilR\n");
        exit(EXIT_FAILURE);
    }

    f->fraction_uptake = malloc(p->core * sizeof(double));
    if (f->fraction_uptake == NULL) {
        fprintf(stderr, "malloc failed allocating soilR\n");
        exit(EXIT_FAILURE);
    }

    f->ppt_gain = malloc(p->core * sizeof(double));
    if (f->ppt_gain == NULL) {
        fprintf(stderr, "malloc failed allocating ppt_gain\n");
        exit(EXIT_FAILURE);
    }

    f->water_loss = malloc(p->core * sizeof(double));
    if (f->water_loss == NULL) {
        fprintf(stderr, "malloc failed allocating water_loss\n");
        exit(EXIT_FAILURE);
    }

    f->water_gain = malloc(p->core * sizeof(double));
    if (f->water_gain == NULL) {
        fprintf(stderr, "malloc failed allocating water_gain\n");
        exit(EXIT_FAILURE);
    }

    /* Depth to bottom of wet soil layers (m) */
    s->water_frac = malloc(p->core * sizeof(double));
    if (s->water_frac == NULL) {
        fprintf(stderr, "malloc failed allocating water_frac\n");
        exit(EXIT_FAILURE);
    }

    /* Depth to bottom of wet soil layers (m) */
    s->wetting_bot = malloc(p->wetting * sizeof(double));
    if (s->wetting_bot == NULL) {
        fprintf(stderr, "malloc failed allocating wetting_bot\n");
        exit(EXIT_FAILURE);
    }

    /* Depth to top of wet soil layers (m) */
    s->wetting_top = malloc(p->wetting * sizeof(double));
    if (s->wetting_top == NULL) {
        fprintf(stderr, "malloc failed allocating wetting_top\n");
        exit(EXIT_FAILURE);
    }

    /* Depth to top of wet soil layers (m) */
    s->wetting_top = malloc(p->wetting * sizeof(double));
    if (s->wetting_top == NULL) {
        fprintf(stderr, "malloc failed allocating wetting_top\n");
        exit(EXIT_FAILURE);
    }

    f->est_evap = malloc(p->core * sizeof(double));
    if (f->est_evap == NULL) {
        fprintf(stderr, "malloc failed allocating est_evap\n");
        exit(EXIT_FAILURE);
    }

    return;
}

void sum_hourly_water_fluxes(fluxes *f, double soil_evap_hlf_hr,
                             double transpiration_hlf_hr, double et_hlf_hr,
                             double interception_hlf_hr,
                             double thoughfall_hlf_hr,
                             double canopy_evap_hlf_hr,
                             double runoff_hlf_hr,
                             double omega_hlf_hr,
                             double rain_hlf_hr) {

    /* add half hour fluxes to day total store */
    f->soil_evap += soil_evap_hlf_hr;
    f->transpiration += transpiration_hlf_hr;
    f->et += et_hlf_hr;
    f->interception += interception_hlf_hr;
    f->throughfall += thoughfall_hlf_hr;
    f->canopy_evap += canopy_evap_hlf_hr;
    f->runoff += runoff_hlf_hr;
    f->omega += omega_hlf_hr; /* average at the end of hour loop */
    f->day_ppt += rain_hlf_hr;

    return;
}

void calc_saxton_stuff(params *p, double *fsoil) {
    /*
        Calculate the key parameters of the Saxton equations:
        cond1, 2, 3 and potA, B

        NB: Currently I'm assuming a single texture for the entire root zone,
            we could set it by layer obviously and we probably should...

        Reference:
        ---------
        * Saxton, K. E., Rawls, W. J., Romberger, J. S. & Papendick, R. I.
          (1986) Estimating generalized soil-water characteristics from texture.
          Soil Science Society of America Journal, 90, 1031-1036.
    */
    int    i;
    double mult1 = 100.0;
    double mult2 = 2.778E-6;
    double mult3 = 1000.0;
    double A = -4.396;
    double B = -0.0715;
    double CC = -4.880E-4;
    double D = -4.285E-5;
    double E = -3.140;
    double F = -2.22E-3;
    double G = -3.484E-5;
    double H = 0.332;
    double J = -7.251E-4;
    double K = 0.1276;
    double P = 12.012;
    double Q = -7.551E-2;
    double R = -3.895;
    double T = 3.671e-2;
    double U = -0.1103;
    double V = 8.7546E-4;
    double x1 = 0.1;
    double x2 = 0.7;
    double tol = 0.0001;
    double dummy = 0.0;
    double sand = fsoil[SAND] * 100.0;
    double clay = fsoil[CLAY] * 100.0;

    /*
    ** As we aren't currently changing texture by layer, the loop is redundant,
    ** but it is probably best to leave it under the assumption we change
    ** this later
    */

    for (i = 0; i < p->core; i++) {
        p->potA[i] = exp(A + B * clay + CC * sand * \
                         sand + D * sand * sand * \
                         clay) * 100.0;
        p->potB[i]  = E + F * clay * clay + G * sand * sand * clay;
        p->cond1[i] = mult2;
        p->cond2[i] = P + Q * sand;
        p->cond3[i] = R + T * sand + U * clay + V * clay * clay;

        p->porosity[i] = H + J * sand + K * log10(clay);

        // field capacity is water content at which SWP = -10 kPa
        //p->field_capacity[i] = zbrent(&saxton_field_capacity, x1, x2, tol,
        //                              p->potA[i], p->potB[i],
        //                              dummy, dummy, dummy);

        // Not sure why Mat used zbrent above, this gives basically the same
        // answer - Saxton '86
        p->field_capacity[i] = exp((2.302 - log(p->potA[i])) / p->potB[i]);
    }

    return;
}

double saxton_field_capacity(double xval, double potA, double potB,
                             double dummy1, double dummy2, double dummy3) {
    /* field capacity calculations for saxton eqns */

    double air_entry_swp = 10.0; /* (kPa) */
    double MPa_2_kPa = 1000.0;
    double swp;

    /* calculate the soil water potential (MPa) */
    swp = -0.001 * potA * pow(xval, potB);
    return (MPa_2_kPa * swp + air_entry_swp);
}

double calc_soil_conductivity(double water_frac, double cond1, double cond2,
                              double cond3) {
    //
    // Soil hydraulic conductivity (m s-1 ) per soil layer based on
    // Saxton et al. (1986) equations. Used in the soil drainage integrator
    double scond;

    /* Avoid floating-underflow error */
    if (water_frac < 0.05) {
        scond = 1E-30;
    } else {
        scond = cond1 * exp(cond2 + cond3 / water_frac);
    }
    return (scond);
}

void calc_soil_water_potential(fluxes *f, params *p, state *s) {
    //
    // Calculate the SWP (MPa) in each soil layer based on algorithms from
    // Saxton et al. (1986). We are not updating the water fraction in each
    // layer
    //

    int    i;
    double arg1, arg2;

    for (i = 0; i < s->rooted_layers; i++) {

        if (s->water_frac[i] > 0.0) {
            f->swp[i] = -0.001 * p->potA[i] * pow(s->water_frac[i], p->potB[i]);
        } else {
            f->swp[i] = -9999.0;
        }
    }

    return;
}

void calc_soil_root_resistance(fluxes *f, params *p, state *s) {

    /* head of pressure (MPa/m) */
    double head = 0.009807;
    double Lsoil, rs, soilR1, soilR2, arg1, arg2, rsum;
    double root_xsec_area = M_PI * p->root_radius * p->root_radius;
    int    i;

    // Store each layers resistance, used in LWP calculatons
    rsum = 0.0;
    for (i = 0; i < s->rooted_layers; i++) {

        /* converts from ms-1 to m2 s-1 MPa-1 */
        Lsoil = f->soil_conduct[i] / head;

        if (Lsoil < 1e-35) {
            /* prevent floating point error */
            f->soilR[i] = 1e35;
        } else {
            rs = sqrt(1.0 / (s->root_length[i] * M_PI));
            arg1 = log(rs / p->root_radius);
            arg2 = 2.0 * M_PI * s->root_length[i] * s->thickness[i] * Lsoil;

            /* soil resistance, convert from MPa s m2 m-3 to MPa s m2 mmol-1 */
            soilR1 = (arg1 / arg2) * 1E-6 * 18. * 0.001;

            // Need to combine resistances in parallel, but we only want the
            // soil term as the root component is part of the plant resistance
            rsum += 1.0 / soilR1;

            // second component of below ground resistance related to root
            // hydraulics
            soilR2 = p->root_resist / (s->root_mass[i] * s->thickness[i]);
            //f->soilR[i] = soilR1 + soilR2; /* MPa s m2 mmol-1 */
            f->soilR[i] = soilR1 + soilR2; /* MPa s m2 mmol-1 */
        }
    }

    f->total_soil_resist = 1.0 / rsum;

    return;
}


void calc_water_uptake_per_layer(fluxes *f, params *p, state *s) {
    //
    // Determine which layer water is extracted from. This is achieved by
    // roughly estimating the maximum rate of water supply from each rooted
    // soil layer, using SWP and hydraulic resistance of each layer. Actual
    // water from each layer is determined using the estimated value as a
    // weighted factor.
    //

    int    i;
    double total_est_evap;

    total_est_evap = 0.0;
    s->weighted_swp = 0.0;

    // Estimate max transpiration from gradient-gravity / soil resistance
    for (i = 0; i < s->rooted_layers; i++) {
        f->est_evap[i] = MAX(0.0, (f->swp[i] - p->min_lwp) / f->soilR[i]);
        total_est_evap += f->est_evap[i];
    }

    if (total_est_evap > 0.0) {
        /* fraction of water taken from layer */
        for (i = 0; i < s->rooted_layers; i++) {
            s->weighted_swp += f->swp[i] * f->est_evap[i];
            f->fraction_uptake[i] = f->est_evap[i] / total_est_evap;
        }
        s->weighted_swp /= total_est_evap;
    } else {
        /* No water was evaporated */
        for (i = 0; i < s->rooted_layers; i++) {
            f->fraction_uptake[i] = 1.0 / (double)s->rooted_layers;
        }
    }

    if (f->fraction_uptake[0] > 1 || f->fraction_uptake[0] < 0) {
        fprintf(stderr, "Problem with the uptake fraction\n");
        exit(EXIT_FAILURE);
    }

    return;
}

void calc_wetting_layers(fluxes *f, params *p, state *s, double soil_evap,
                         double surface_water) {
    //
    // Tracks surface wetting and drying in the top soil layer and so the
    // thickness of the uppermost dry layer and thus soil evaporation
    //
    // wetting_bot - Depth to bottom of wet soil layers (m)
    // wetting_top - Depth to top of wet soil layers (m)
    //

    double seconds_per_step = 1800.0;
    double dmin = 0.001;
    double airspace = p->porosity[0];
    double min_val, netc, diff;
    int    i, ar1, ar2;

    //
    // soil LE should be withdrawn from the wetting layer with the
    // smallest depth..
    //
    ar1 = 0;
    min_val = 9999.9;
    for (i = 0; i < p->wetting; i++) {
        if (s->wetting_bot[i] > 0.0 && s->wetting_bot[i] < min_val) {
            ar1 = i;
            min_val = s->wetting_bot[i];
        }
    }

    // Need to make this into a negative flux of energy from the soil's
    // perspective to match remaining logic here, i.e. negative leaving the
    // surface
    if (soil_evap > 0.0)
        soil_evap *= -1.0;

    // Calulate the net change in wetting in the top zone
    netc = (soil_evap * MM_TO_M) / airspace + \
           (surface_water * MM_TO_M) / airspace;

    // wetting
    if (netc > 0.0) {

        /*
        ** resaturate the layer if top is dry and recharge is greater
        **  than dry_thick
        */

        if ((netc > s->wetting_top[ar1]) && (s->wetting_top[ar1] > 0.0)) {

            // extra water to deepen wetting layer
            diff = netc - s->wetting_top[ar1];
            s->wetting_top[ar1] = 0.0;
            if (ar1 > 0) {
                // Not in primary layer (primary layer can't extend deeper)
                s->wetting_bot[ar1] += diff;
            }
            s->dry_thick = dmin;
        } else {

            if (s->wetting_top[ar1] == 0.0) {

                // surface is already wet, so extend depth of this wet zone
                if (ar1 > 0) {
                    // not in primary lay (primary layer can't extend deeper)
                    s->wetting_bot[ar1] += netc;
                    if (s->wetting_bot[ar1] >= s->wetting_top[ar1-1]) {
                        // Layers are conterminous..
                        s->wetting_top[ar1-1] = s->wetting_top[ar1];
                        s->wetting_top[ar1] = 0.;     /* remove layer */
                        s->wetting_bot[ar1] = 0.;    /* remove layer */
                    }
                }
            } else {

                // Create a new wetting zone
                s->wetting_top[ar1+1] = 0.0;
                s->wetting_bot[ar1+1] = netc;
            }
            s->dry_thick = dmin;
        }

    // Drying
    } else {
        // Drying increases the depth to top of wet soil layers
        s->wetting_top[ar1] -= netc;

        // Wetting layer is dried out.
        if (s->wetting_top[ar1] >= s->wetting_bot[ar1]) {
            /* How much more drying is there? */
            diff = s->wetting_top[ar1] - s->wetting_bot[ar1];
            s->wetting_top[ar1] = 0.0;
            s->wetting_bot[ar1] = 0.0;
            ar2 = (ar1 + 1) - 1;          // +1 is fortran logic offset

            /* Move to deeper wetting layer */
            if (ar2 > 0) {
                /* dry out deeper layer */
                s->wetting_top[ar2] += diff;
                s->dry_thick = MAX(dmin, s->wetting_top[ar2]);
            /* no deeper layer */
            } else {
                /* layer 1 is dry */
                s->dry_thick = s->thickness[0];
            }
        } else {
            s->dry_thick = MAX(dmin, s->wetting_top[ar1]);
        }
    }

    if (s->dry_thick == 0.0) {
        fprintf(stderr, "Problem in dry_thick\n");
        exit(EXIT_FAILURE);
    }

    return;
}

double calc_infiltration(fluxes *f, params *p, state *s, double surface_water) {
    /*
        Takes surface_water and distrubutes it among top layers. Assumes
        total infilatration in timestep.
    */
    int    i;
    double add;     // surface water available for infiltration (m)
    double wdiff;   // available space in a given soil lyr for water to fill (m)
    double runoff;

    add = surface_water * MM_TO_M;

    for (i = 0; i < p->core; i++) {
        f->ppt_gain[i] = 0.0;
    }

    runoff = 0.0;
    for (i = 0; i < p->soil_layers; i++) {

        // determine the available pore space in current soil layer
        wdiff = MAX(0.0, (p->porosity[i] - s->water_frac[i]) * \
                          s->thickness[i] - f->water_gain[i] + \
                          f->water_loss[i]);

        // is the input of water greater than available space?
        if (add > wdiff) {
            // if so fill and subtract from input and move on to the next layer
            f->ppt_gain[i] = wdiff;
            add -= wdiff;
        } else {
            // otherwise infiltate all in the current layer
            f->ppt_gain[i] = add;
            add = 0.0;
        }

        // if we have added all available water we are done
        if (add <= 0.0) {
            break;
        }
    }

    // if after all of this we have some water left assume it is runoff
    if (add >  0.0) {
       runoff = add * M_TO_MM;
    } else {
       runoff = 0.0;
    }


    return (runoff);
}


void calc_soil_balance(fluxes *f, nrutil *nr, params *p, state *s,
                       int soil_layer, double *water_lost) {
    //
    // Integrator for soil gravitational drainage
    //
    int    nbad;                /* N of unsuccessful changes of the step size */
    int    nok;                 /* N of successful changes of the step size */
    int    i, N = 1, max_iter;
    double eps = 1.0e-4;        /* precision */
    double h1 = .001;           /* first guess at integrator size */
    double hmin = 0.0;          /* minimum value of the integrator step */
    double x1 = 1.0;             /* initial time */
    double x2 = 2.0;             /* final time */

    /* value affecting the max time interval at which variables should b calc */
    double  soilpor = p->porosity[soil_layer];
    double  unsat, drain_layer, new_water_frac, change;
    //double *ystart = NULL;
    //ystart = dvector(1,N);

    for (i = 1; i <= N; i++) {
        nr->ystart[i] = 0.0;
    }

    /* unsaturated volume of layer below (m3 m-2) */
    unsat = MAX(0.0, (p->porosity[soil_layer+1] - \
                      s->water_frac[soil_layer+1]) *\
                      s->thickness[soil_layer+1] / s->thickness[soil_layer]);

    /* soil water capacity of the current layer */
    drain_layer = p->field_capacity[soil_layer];

    // SPA assumes drainages occurs if water_frac > field_capacity, but we are
    // assuming drainage always takes place. Remko argues for this making more
    // physical sense
    if (s->water_frac[soil_layer] > 0.0) {

        // ystart is a vector 1..N, so need to index from 1 not 0
        nr->ystart[1] = s->water_frac[soil_layer];

        // Runge-Kunte ODE integrator used to estimate soil gravitational
        // drainage during each time-step
        odeint(nr->ystart, N, x1, x2, eps, h1, hmin, &nok, &nbad, unsat,
               drain_layer, p->cond1[soil_layer], p->cond2[soil_layer],
               p->cond3[soil_layer], nr, soil_water_store, rkqs);

        /* ystart is a vector 1..N, so need to index from 1 */
        new_water_frac = nr->ystart[1];

        /* convert from water fraction to absolute amount (m) */
        change = (s->water_frac[soil_layer] - new_water_frac) * \
                  s->thickness[soil_layer];

        // update soil layer below with drained liquid
        f->water_gain[soil_layer+1] += change;
        f->water_loss[soil_layer] += change;

        // SPA assumption. Water can be passed through the final layer
        // (to the core layer), but this water is lost in any balance as this
        // layer is actually kept dry to ensure drainage occurs.
        // I'm going to add this water to runoff just to ensure water balances.
        // I guess it can be thought of as sub-surface runoff

        // Counting from zero so this logic makes sense
        if (soil_layer+1 == p->core - 1) {
            *water_lost += change;
        }
    }

    if (f->water_loss[soil_layer] < 0.0) {
        fprintf(stderr, "waterloss probem in soil_balance: %d %f\n",
                soil_layer, f->water_loss[soil_layer]);
        exit(EXIT_FAILURE);
    }

    //free_dvector(ystart, 1, N);

    return;
}


void soil_water_store(double time_dummy, double y[], double dydt[],
                      double unsat, double drain_layer, double cond1,
                      double cond2, double cond3) {

    /*
    ** numerical lib vectors are index 1..n, so we need to index the return
    ** from the odeint func with 1, not 0
    */
    int index = 1;

    double drainage;

    drainage = calc_soil_conductivity(y[index], cond1, cond2, cond3);

    // Convert units, soil conductivity is in m s-1 //
    drainage *= SEC_2_HLFHR;

    // gravitational drainage above field_capacity
    if (y[index] <= drain_layer) {
        drainage = 0.0;
    }

    // layer below cannot accept more water than unsat
    if (drainage > unsat) {
        drainage = unsat;
    }

    // waterloss from this layer
    dydt[index] = -drainage;

    return;
}

void extract_water_from_layers(fluxes *f, state *s, double soil_evap,
                               double transpiration) {

    // Extract soil evaporation and transpiration from the soil profile

    int rr, i;

    // Is soil evap taken from first or second layer?
    if (s->dry_thick < s->thickness[0]) {
        // The dry zone does not extend beneath the top layer
        rr = 0;
    } else {
        // The dry zone does extend beneath the top layer
        rr = 1;
    }

    if (soil_evap > 0.0) {
      f->water_loss[rr] += soil_evap * MM_TO_M;
    } // ignoring water gain due to due formation...


    // Determing water loss from each layer due to transpiration
    for (i = 0; i < s->rooted_layers; i++) {
        f->water_loss[i] += (transpiration * MM_TO_M) * \
                             f->fraction_uptake[i];
    }

    return;
}

void update_soil_water_storage(fluxes *f, params *p, state *s,
                               double *soil_evap, double *transpiration) {
    //
    // Update the soil water storage at the end of the timestep
    //
    double root_zone_total, water_content, needed, taken, prev_soil_evap;
    int    i, rr;
    double soil_evap_overshoot, transpiration_overshoot, prev_trans;
    double effective_swp, wp;



    root_zone_total = 0.0;
    for (i = 0; i < p->soil_layers; i++) {

        // water content of soil layer (m)
        water_content = s->water_frac[i] * s->thickness[i];

        needed = water_content + f->water_gain[i] + \
                 f->ppt_gain[i] - f->water_loss[i];

        // Is soil evap taken from first or second layer?
        if (s->dry_thick < s->thickness[0]) {
            // The dry zone does not extend beneath the top layer
            rr = 0;
        } else {
            // The dry zone does extend beneath the top layer
            rr = 1;
        }

        // Correction for potential to over-evaporate if using Emax drought
        // stress correction. This stops that happening.
        if (i == rr) {
            if (needed < 0.0) {
                prev_soil_evap = *soil_evap;
                *soil_evap = MAX(0.0, *soil_evap + (needed * M_TO_MM));
                if (*soil_evap > 0.0) {
                    taken = (prev_soil_evap - *soil_evap) * MM_TO_M;
                    needed -= taken;
                    f->water_loss[i] += taken;
                }
            }
        } else {
            if (needed < 0.0) {
                prev_trans = *transpiration;
                *transpiration = MAX(0.0, *transpiration + (needed * M_TO_MM));
                taken = (prev_trans - *transpiration) * MM_TO_M;
                f->water_loss[i] += taken;
            }
        }

        // NB water gain here is drainage from the layer above
        water_content = MAX(0.0, water_content +    \
                                 f->water_gain[i] + \
                                 f->ppt_gain[i] -   \
                                 f->water_loss[i]);

        // Determine volumetric water content water content of layer (m3 m-3)
        s->water_frac[i] = water_content / s->thickness[i];

        // update old GDAY effective two-layer buckets
        // - this is just for outputting, these aren't used.
        if (i == 0) {
            s->pawater_topsoil = water_content * M_TO_MM;

        } else {

            // SPA doesn't have a wilting point per layer, we can infer an
            // effective one by assuming a SWP of -1.5 MPa. I think by doing
            // this our water balance check would no longer work...so expect
            // some error there
            effective_swp = -1.5;

            // I've just rearranged the SWP calculation
            wp = pow( effective_swp / (-0.001 * p->potA[i]),
                      (1.0 / p->potB[i]) );
            root_zone_total += MAX(0.0, (s->water_frac[i] - wp) * \
                                         s->thickness[i] * M_TO_MM);
        }
    }
    s->pawater_root = root_zone_total;



    return;
}

double calc_xylem_water_potential(double rwc, double capac) {
    //
    // Calculate the stem xylem water potential (P), based on relative water
    // content (RWC) and capacitance.
    //
    // Parameters:
    // ----------
    // rwc : double
    //    relative water content [-]
    // capac : double
    //    capacitance (MPa per unit relative water content)
    //
    // Returns:
    // --------
    // xylem_psi : double
    //  xylem water potential (MPa)
    //
    double psi1, psi2, xylem_psi, arg1, arg2, arg3;
    double break0 = 0.5;    // determines shape of asymptote function
    double hmshape = 0.99;  // determines shape of hyperbolic minimum

    // safety
    if (rwc > 1.0) {
        rwc = 1.0;
    }

    // linear dependence over most of the range.
    psi1 = -(1.0 - rwc) / capac;

    // when approaching zero rwc, the water potential has to go to infinity.
    psi2 = -log(rwc / break0);
    if (psi2 < 0.0) {
        psi2 = 0.0;
    }
    psi2 = -psi2;

    // hyperbolic minimum
    arg1 = psi1 + psi2;
    arg2 = (psi1 + psi2) * (psi1 + psi2);
    arg3 = 4.0 * hmshape * psi1 * psi2;
    xylem_psi = (arg1 - sqrt(arg2 - arg3)) / (2.0 * hmshape);

    return (xylem_psi);
}

double calc_relative_weibull(double p, double p50, double sx) {
    //
    // Calculate the relative conductivity, given xylem water potential (p),
    // the p50, and the shape parameter (sx)
    //
    // Parameters:
    // ----------
    // p : double
    //    xylem water potential (MPa)
    // p50 : double
    //    xylem water potential where 50% of the conductivity is lost
    // sx : double
    //    slope paramater: derivative (% MPa-1) at x (e.g. s50 is the slope
    //    of the curve at P50). Higher values thus indicate steeper response
    //    to xylem pressure
    //
    // Reference:
    // ---------
    // * Ogle et al. (2009) Ecological Applications, 19, 577-581.
    //
    double v, relative_weibull;

    v = -50.0 * log(0.5);
    relative_weibull = 1.0 - pow(0.5, pow((p / p50), (p50 * sx) / v));

    return (relative_weibull);
}

void update_plant_water_store(canopy_wk *cw, params *p, state *s,
                              double *transpiration, double *et,
                              double et_deficit, double year, double doy) {

    // 5 % of full hydration
    double min_value = 0.05 * cw->plant_water0;
    double ratio, water_flux, stem_relk, arg1, arg2;
    double delta_water_store = 0.0;
    double conv;

    // Under normal circumstances, i.e et_deficit = 0, the assumption is that
    // they will fill up the stem immediately (even if near empty).
    // stem water potential is soilwp - transpiration / (2*k)
    if (et_deficit * MOLE_WATER_2_G_WATER * SEC_2_HLFHR < 1E-06) {

        // mm 30 min-1 -> mmol m-2 s-1
        conv = KG_AS_G * G_WATER_2_MOL_WATER * MOL_2_MMOL * HLFHR_2_SEC;
        water_flux = *transpiration * conv;

        arg1 = s->weighted_swp;
        arg2 = water_flux / (2.0 * cw->plant_k * s->lai);
        cw->xylem_psi = arg1 - arg2;

        // refill plant water store
        cw->plant_water = cw->plant_water0 * (1.0 + cw->xylem_psi * p->capac);

    } else {

        // now reduce stem water content even further by amount of
        // transpiration that is not sustained by soil water uptake
        conv = MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
        cw->plant_water -= et_deficit * conv;

        // To avoid stopping the simulation when we are "dead"
        // (i.e. xylempsi is very low)
        if (cw->plant_water < min_value) {
            cw->plant_water = min_value;
        }

        // and recalculate corresponding xylem water potential
        ratio = cw->plant_water / cw->plant_water0;
        cw->xylem_psi = calc_xylem_water_potential(ratio, p->capac);

    }

    // Need to add water we took from the plant store to transpiration output
    //
    // mol m-2 s-1 to mm/30min
    conv = MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
    *transpiration += et_deficit * conv;
    *et += et_deficit * conv;

    // stem relative conductivity (0-1)
    stem_relk = calc_relative_weibull(cw->xylem_psi, p->p50, p->plc_shape);

    // is the plant dead? Going to store this information and we can work
    // out how to write it to a file later.
    if (stem_relk < (1.0 - p->plc_dead) && cw->not_dead) {
        cw->not_dead = FALSE;
        cw->death_year = (double)year;
        cw->death_doy = (double)doy;
    }

    // Update plant conductance
    cw->plant_k = stem_relk * p->kp;

    return;
}

void calc_soil_balance_cascading(fluxes *f, nrutil *nr, params *p, state *s,
                                 int soil_layer, double *water_lost) {
    //
    // Much simpler solution to soil water drainge: assumes that water can move
    // (downwards) through the soil profile, filling up the layers until
    // field capacity is reached, with the fraction of water exceeding this
    // threshold moving to the deeper layer
    //

    int     i;
    double  unsat, drain_layer, liquid, new_water_frac, change, drainage;

    /* unsaturated volume of layer below (m3 m-2) */
    unsat = MAX(0.0, (p->porosity[soil_layer+1] - \
                      s->water_frac[soil_layer+1]) *\
                      s->thickness[soil_layer+1] / s->thickness[soil_layer]);

    /* soil water capacity of the current layer */
    drain_layer = p->field_capacity[soil_layer];
    liquid = s->water_frac[soil_layer];

    if (liquid > 0.0 && liquid > drain_layer) {

        // Assumption that ~5% of the water in the layer will drain, this is
        // taken from running the standard SPA model and averaging the change
        // for 10000 timesteps. I got ~7% for the site in question, tumba, so
        // I've assumed a universal 10%, this number could do with more
        // testing :)

        // waterloss from this layer
        drainage = liquid * 0.1;

        // gravitational drainage above field_capacity
        if (drainage <= drain_layer) {
            drainage = 0.0;
        }

        // layer below cannot accept more water than unsat
        if (drainage > unsat) {
            drainage = unsat;
        }

        // waterloss from this layer
        new_water_frac = s->water_frac[soil_layer] - drainage;

        /* convert from water fraction to absolute amount (m) */
        change = (s->water_frac[soil_layer] - new_water_frac) * \
                  s->thickness[soil_layer];

        // Update current layer
        f->water_loss[soil_layer] += change;

        // update soil layer below with drained liquid
        if (soil_layer+1 < p->soil_layers) {
            f->water_gain[soil_layer+1] += change;
        } else {
            // We are draining through the bottom soil layer, add to runoff
            *water_lost += change;
        }

    }

    return;
}
