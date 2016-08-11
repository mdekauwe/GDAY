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
        setup_hydraulics_arrays(f, p, s);

        for (i = 0; i < p->wetting; i++) {
            s->wetting_bot[i] = 0.0;
            s->wetting_top[i] = 0.0;
        }

        /* saturate the top layer */
        s->wetting_bot[0] = s->thickness[0];

        /* Initalise SW fraction - we should read this from param file */
        s->initial_water = 0.0;
        for (i = 0; i < p->n_layers; i++) {
            s->water_frac[i] = 0.4;
            s->initial_water += 1E3 * (s->water_frac[i] * s->thickness[i]);
        }

    }

    free(fsoil_top);
    free(fsoil_root);

    return;
}

void calculate_water_balance_sub_daily(control *c, fluxes *f, met *m,
                                       params *p, state *s, int daylen,
                                       double trans_leaf, double omega_leaf,
                                       double rnet_leaf) {
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
        length of day in hours. (Dummy argument, only passed for daily model)
    trans_leaf : double
        leaf transpiration (Dummy argument, only passed for sub-daily model)
    omega_leaf : double
        decoupling coefficient (Dummy argument, only passed for sub-daily model)
    rnet_leaf : double
        total canopy rnet (Dummy argument, only passed for sub-daily model)

    */
    int    i, rr;
    double soil_evap, et, interception, runoff, conv, transpiration, net_rad;
    double SEC_2_DAY, DAY_2_SEC, transpiration_am, transpiration_pm, gs_am;
    double gs_pm, LE_am, LE_pm, ga_am, ga_pm, net_rad_am, net_rad_pm, omega_am;
    double gpp_am, gpp_pm, omega_pm, throughfall, canopy_evap;
    double water_content, surface_water;



    if (c->water_balance == HYDRAULICS) {

        /* assume that water infiltrates within a timestep */
        surface_water = 0.0;

        /*
        ** The loop needs to be outside the func as we need to be able to
        ** calculate the soil conductance per layer and call this via
        ** the integration func when we update the soil water balance
        */
        for (i = 0; i < p->n_layers; i++) {
            f->soil_conduct[i] = calc_soil_conductivity(s->water_frac[i],
                                                        p->cond1[i], p->cond2[i],
                                                        p->cond3[i]);
        }
        calc_soil_water_potential(f, p, s);
        calc_soil_root_resistance(f, p, s);
        calc_water_uptake_per_layer(f, p, s);


        /* calculate potential canopy evap rate, this may be reduced later
           depending on canopy water storage */
        canopy_evap = calc_canopy_evaporation(m, p, s, rnet_leaf);

        /* mol m-2 s-1 to mm d-1 */
        conv = MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
        canopy_evap *= conv;

        /* We could now replace this interception bit with the Rutter scheme? */
        calc_interception(c, m, p, f, s, &throughfall, &interception,
                          &canopy_evap);

        net_rad = calc_net_radiation(p, m->sw_rad, m->tair);
        soil_evap = calc_soil_evaporation(m, p, s, net_rad);
        soil_evap *= MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;

        /* mol m-2 s-1 to mm/30 min */
        transpiration = trans_leaf * MOLE_WATER_2_G_WATER * G_TO_KG * \
                        SEC_2_HLFHR;

        et = transpiration + soil_evap + canopy_evap;

        /* update surface water puddle */
        surface_water += throughfall;

        /* determine water loss in upper layers due to evaporation */
        calc_wetting_layers(f, p, s, soil_evap, surface_water);

        /*  From which layer is evap water withdrawn? */
        if (s->dry_thick < s->thickness[0]) {
            rr = 1;     /* The dry zone does not extend beneath the top layer */
        } else {
            rr = 2;     /* The dry zone does extend beneath the top layer */
        }

        /* determine water loss in upper layers due to evaporation */
        if (soil_evap > 0.0) {
          /* Evaporation (m 30min-1) */
          f->water_loss[rr] += soil_evap * MM_TO_M;
        } /* ignoring water gain due to due formation...


        /* Determing water loss from each layer due to transpiration */
        for (i = 0; i < s->rooted_layers; i++) {
            f->water_loss[i] += (transpiration * MM_TO_M) * \
                                    f->fraction_uptake[i];
        }

        /*
        ** determines water movement between soil layers due drainage
        ** down the profile
        */
        for (i = 0; i < p->n_layers; i++) {
            calc_soil_balance(f, p, s, i);
        }

        /*
        ** how much surface water infiltrantes the first soil layer in the
        ** current time step? Water which does not infiltrate in a single step
        ** is considered runoff
        */
        runoff = calc_infiltration(f, p, s, surface_water);
        calc_soil_water_potential(f, p, s);
        calc_soil_root_resistance(f, p, s);

        /* Update the soil water storage */
        for (i = 0; i < p->n_layers; i++) {

            /* water content of soil layer (m) */
            water_content = s->water_frac[i] * s->thickness[i];

            /*
            ** Net change in water content (m);
            ** max condition to ensure
            */
            water_content = MAX(0.0, water_content + \
                                     f->water_gain[i] + \
                                     f->ppt_gain[i] - \
                                     f->water_loss[i]);

            /* Determine new total water content of layer (m) */
            s->water_frac[i] = water_content / s->thickness[i];

        }
        exit(1);


    } else {

        /* Simple soil water bucket appoximation */

        /* calculate potential canopy evap rate, this may be reduced later
           depending on canopy water storage */
        canopy_evap = calc_canopy_evaporation(m, p, s, rnet_leaf);

        /* mol m-2 s-1 to mm/day */
        conv = MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
        canopy_evap *= conv;
        calc_interception(c, m, p, f, s, &throughfall, &interception,
                          &canopy_evap);

        net_rad = calc_net_radiation(p, m->sw_rad, m->tair);
        soil_evap = calc_soil_evaporation(m, p, s, net_rad);
        soil_evap *= MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;

        /* mol m-2 s-1 to mm/30 min */
        transpiration = trans_leaf * MOLE_WATER_2_G_WATER * G_TO_KG * \
                        SEC_2_HLFHR;

        /*
        ** NB. et, transpiration & soil evap may all be adjusted in
        ** update_water_storage if we don't have sufficient water
        */
        et = transpiration + soil_evap + canopy_evap;

        update_water_storage(c, f, p, s, throughfall, interception, canopy_evap,
                             &transpiration, &soil_evap, &et, &runoff);

    }

    sum_hourly_water_fluxes(f, soil_evap, transpiration, et, interception,
                            throughfall, canopy_evap, runoff, omega_leaf);
}


void setup_hydraulics_arrays(fluxes *f, params *p, state *s) {
    /* Allocate the necessary memory for all the hydraulics arrays */

    f->soil_conduct = malloc(p->n_layers * sizeof(double));
    if (f->soil_conduct == NULL) {
        fprintf(stderr, "malloc failed allocating soil_conduct\n");
        exit(EXIT_FAILURE);
    }

    f->swp = malloc(p->n_layers * sizeof(double));
    if (f->swp == NULL) {
        fprintf(stderr, "malloc failed allocating swp\n");
        exit(EXIT_FAILURE);
    }

    f->soilR = malloc(p->n_layers * sizeof(double));
    if (f->soilR == NULL) {
        fprintf(stderr, "malloc failed allocating soilR\n");
        exit(EXIT_FAILURE);
    }

    f->fraction_uptake = malloc(p->n_layers * sizeof(double));
    if (f->fraction_uptake == NULL) {
        fprintf(stderr, "malloc failed allocating soilR\n");
        exit(EXIT_FAILURE);
    }

    f->ppt_gain = malloc(p->n_layers * sizeof(double));
    if (f->ppt_gain == NULL) {
        fprintf(stderr, "malloc failed allocating ppt_gain\n");
        exit(EXIT_FAILURE);
    }

    f->water_loss = malloc(p->n_layers * sizeof(double));
    if (f->water_loss == NULL) {
        fprintf(stderr, "malloc failed allocating water_loss\n");
        exit(EXIT_FAILURE);
    }

    f->water_gain = malloc(p->n_layers * sizeof(double));
    if (f->water_gain == NULL) {
        fprintf(stderr, "malloc failed allocating water_gain\n");
        exit(EXIT_FAILURE);
    }


    /* Depth to bottom of wet soil layers (m) */
    s->water_frac = malloc(p->n_layers * sizeof(double));
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

}



void sum_hourly_water_fluxes(fluxes *f, double soil_evap_hlf_hr,
                             double transpiration_hlf_hr, double et_hlf_hr,
                             double interception_hlf_hr,
                             double thoughfall_hlf_hr,
                             double canopy_evap_hlf_hr,
                             double runoff_hlf_hr, double omega_hlf_hr) {

    /* add half hour fluxes to day total store */
    f->soil_evap += soil_evap_hlf_hr;
    f->transpiration += transpiration_hlf_hr;
    f->et += et_hlf_hr;
    f->interception += interception_hlf_hr;
    f->throughfall += thoughfall_hlf_hr;
    f->canopy_evap += canopy_evap_hlf_hr;
    f->runoff += runoff_hlf_hr;
    f->omega += omega_hlf_hr; /* average at the end of hour loop */

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

    p->potA = malloc(p->n_layers * sizeof(double));
    if (p->potA == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's potA\n");
        exit(EXIT_FAILURE);
    }

    p->potB = malloc(p->n_layers * sizeof(double));
    if (p->potB == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's potB\n");
        exit(EXIT_FAILURE);
    }

    p->cond1 = malloc(p->n_layers * sizeof(double));
    if (p->cond1 == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's cond1\n");
        exit(EXIT_FAILURE);
    }

    p->cond2 = malloc(p->n_layers * sizeof(double));
    if (p->cond1 == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's cond2\n");
        exit(EXIT_FAILURE);
    }

    p->cond3 = malloc(p->n_layers * sizeof(double));
    if (p->cond1 == NULL) {
        fprintf(stderr, "malloc failed allocating Saxton's cond3\n");
        exit(EXIT_FAILURE);
    }

    p->porosity = malloc(p->n_layers * sizeof(double));
    if (p->porosity == NULL) {
        fprintf(stderr, "malloc failed allocating porosity\n");
        exit(EXIT_FAILURE);
    }

    p->field_capacity = malloc(p->n_layers * sizeof(double));
    if (p->field_capacity == NULL) {
        fprintf(stderr, "malloc failed allocating field_capacity\n");
        exit(EXIT_FAILURE);
    }

    /*
    ** As we aren't currently changing texture by layer, the loop is redundant,
    ** but it is probably best to leave it under the assumption we change
    ** this later
    */


    for (i = 0; i < p->n_layers; i++) {
        p->potA[i] = exp(A + B * clay + CC * sand * \
                         sand + D * sand * sand * \
                         clay) * 100.0;
        p->potB[i]  = E + F * clay * clay + G * sand * sand * clay;
        p->cond1[i] = mult2;
        p->cond2[i] = P + Q * sand;
        p->cond3[i] = R + T * sand + U * clay + V * clay * clay;

        p->porosity[i] = H + J * sand + K * log10(clay);

        /* field capacity is water content at which SWP = -10 kPa */
        p->field_capacity[i] = zbrent(&saxton_field_capacity, x1, x2, tol,
                                      p->potA[i], p->potB[i],
                                      dummy, dummy, dummy);

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
    /* Soil conductivity (m s-1 ) per layer */
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
    /*
        Calculate the SWP (MPa) without updating the water fraction in each
        layer
    */
    int    i;
    double arg1, arg2;

    for (i = 0; i < p->n_layers; i++) {

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

    /* saxton water retention equation are off by default */
    double ab = 1.0;
    double Lsoil, rs, rs2, soilR1, soilR2;
    int    i;

    double root_xsec_area = M_PI * p->root_radius * p->root_radius;

    for (i = 0; i < p->n_layers; i++) {
        /* converts from ms-1 to m2 s-1 MPa-1 */
        Lsoil = f->soil_conduct[i] / head;

        if (Lsoil < 1e-35) {
            /* prevent floating point error */
            f->soilR[i] = 1e35;
        } else {

            rs = sqrt(1.0 / (s->root_length[i] * M_PI));
            rs2 = log(rs / p->root_radius) / \
                     (2.0 * M_PI * s->root_length[i] * s->thickness[i] * Lsoil);

            /* soil resistance, convert from MPa s m2 m-3 to MPa s m2 mmol-1 */
            soilR1 = rs2 * 1E-6 * 18 * 0.001;

            /*
            ** second component of below ground resistance related to root
            ** hydraulics
            */
            soilR2 = p->root_resist / (s->root_mass[i] * s->thickness[i] / ab);
            f->soilR[i] = soilR1 + soilR2; /* MPa s m2 mmol-1 */
        }
    }

    return;
}


void calc_water_uptake_per_layer(fluxes *f, params *p, state *s) {

    /* Figure out which layer water is extracted from */
    int    i;
    double est_evap[p->n_layers];
    double total_est_evap;

    /* Estimate max transpiration from gradient-gravity / soil resistance. */
    total_est_evap = 0.0;
    for (i = 0; i < s->rooted_layers; i++) {
        est_evap[i] = MAX(0.0, (f->swp[i] - p->min_lwp) / f->soilR[i]);
        total_est_evap += est_evap[i];
    }

    /* Water was evaporated from some layers..*/
    s->weighted_swp = 0.0;
    if (total_est_evap > 0.0) {
        for (i = 0; i < p->n_layers; i++) {
            s->weighted_swp += f->swp[i] * est_evap[i];
            /* fraction of water taken from layer */
            f->fraction_uptake[i] = est_evap[i] / total_est_evap;
        }
        s->weighted_swp /= total_est_evap;
    } else {
        /* No water was evaporated */
        f->fraction_uptake[i] = 1.0 / (double)p->n_layers;
    }

    if (f->fraction_uptake[0] > 1 || f->fraction_uptake[0] < 0) {
        fprintf(stderr, "Problem with the uptake fraction\n");
        exit(EXIT_FAILURE);
    }

    return;
}

void calc_wetting_layers(fluxes *f, params *p, state *s, double soil_evap,
                         double surface_water) {
    /*
        surface wetting and drying determines thickness of dry layer
        and thus soil evaporation
    */

    double seconds_per_step = 1800.0;
    double dmin = 0.001;
    double airspace = p->porosity[0];
    double min_val, netc, diff;
    int    i, ar1, ar2;

    /*
    ** soil LE should be withdrawn from the wetting layer with the
    ** smallest depth..
    */
    ar1 = 0;
    min_val = 9999.9;
    for (i = 0; i < p->n_layers; i++) {
        if (s->wetting_bot[i] > 0.0) {
            if (s->wetting_bot[i] < min_val) {
                ar1 = i;
                min_val = s->wetting_bot[i];
            }
        }
    }

    /* Calulate the net change in wetting in the top zone */
    netc = (soil_evap * MM_TO_M) / airspace + \
           (surface_water * MM_TO_M) / airspace;


    /* wetting */
    if (netc > 0.0) {
        /*
        ** resaturate the layer if top is dry and recharge is greater
        **  than dry_thick
        */
        if ((netc > s->wetting_top[ar1]) && (s->wetting_top[ar1] > 0.0)) {
            /* extra water to deepen wetting layer */
            diff = netc - s->wetting_top[ar1];
            s->wetting_top[ar1] = 0.0;
            if (ar1 > 0) {
                /* Not in primary layer (primary layer can't extend deeper) */
                s->wetting_bot[ar1] += diff;
            }
            s->dry_thick = dmin;
        } else {
            if (s->wetting_top[ar1] == 0.0) {
                /* surface is already wet, so extend depth of this wet zone */
                if (ar1 > 0) {
                    /* not in primary lay (primary layer can't extend deeper) */
                    s->wetting_bot[ar1] += netc;
                    if (s->wetting_bot[ar1] >= s->wetting_top[ar1-1]) {
                        /* Layers are conterminous.. */
                        s->wetting_top[ar1-1] = s->wetting_top[ar1];
                        s->wetting_top[ar1] = 0.;     /* remove layer */
                        s->wetting_bot[ar1] = 0.;    /* remove layer */
                    }
                }
            } else {
                /* Create a new wetting zone */
                s->wetting_top[ar1+1] = 0.0;
                s->wetting_bot[ar1+1] = netc;
            }
            s->dry_thick = dmin;
        }
    /* drying */
    } else {
        /* Drying increases the wettingtop depth */
        s->wetting_top[ar1] -= netc;

        /* Wetting layer is dried out. */
        if (s->wetting_top[ar1] > s->wetting_bot[ar1]) {
            /* How much more drying is there? */
            diff = s->wetting_top[ar1] - s->wetting_bot[ar1];
            s->wetting_top[ar1] = 0.0;
            s->wetting_bot[ar1] = 0.0;
            ar2 = ar1;
            ar2 -= 1;

            /* Move to deeper wetting layer */
            if (ar2 > 0.0) {
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

double calc_infiltration(fluxes *f, params *p, state *s,
                         double surface_water) {
    /*
        Takes surface_water and distrubutes it among top layers. Assumes
        total infilatration in timestep.
    */
    int    i;
    double add, wdiff, runoff;

    add = surface_water * MM_TO_M;

    for (i = 0; i < p->n_layers; i++) {
        f->ppt_gain[i] = 0.0;
    }

    runoff = 0.0;
    for (i = 0; i < p->n_layers; i++) {
        wdiff = MAX(0.0, (p->porosity[i] - s->water_frac[i]) * \
                          s->thickness[i] - f->water_gain[i] + \
                          f->water_loss[i]);
        if (add > wdiff) {
            f->ppt_gain[i] = wdiff;
            add -= wdiff;
        } else {
            f->ppt_gain[i] = add;
            add = 0.0;
        }

        if (add < 0.0) {
            fprintf(stderr, "Error in infiltration calculation\n");
            exit(EXIT_FAILURE);
        }

        if (add > 0.0) {
            runoff += add;
        }
    }

    return (runoff);
}


void calc_soil_balance(fluxes *f, params *p, state *s, int soil_layer) {
    /* Integrator for soil gravitational drainage */

    int    nbad;                /* N of unsuccessful changes of the step size */
    int    nok;                 /* N of successful changes of the step size */
    int    N = 1, max_iter;
    double eps = 1.0e-3;        /* precision */
    double h1 = .001;           /* first guess at integrator size */
    double hmin = 0.0;          /* minimum value of the integrator step */
    double x1 = 1.0;             /* initial time */
    double x2 = 2.0;             /* final time */
    /* value affecting the max time interval at which variables should b calc */
    double soilpor = p->porosity[soil_layer];
    double unsat, drain_layer, liquid, new_water_frac, change;

    double *ystart;
    ystart = dvector(1,N);

    /* unsaturated volume of layer below (m3 m-2) */
    unsat = MAX(0.0, (p->porosity[soil_layer+1] - \
                      s->water_frac[soil_layer+1]) *\
                      s->thickness[soil_layer+1] / s->thickness[soil_layer]);

    /* soil water capacity of the current layer */
    drain_layer = p->field_capacity[soil_layer];
    liquid = s->water_frac[soil_layer];

    /*
    ** initial conditions; i.e. is there liquid water and more
    ** water than layer can hold
    */
    if ( (liquid > 0.0) && (s->water_frac[soil_layer] > drain_layer) ) {

        /* there is liquid water */

        /* ystart is a vector 1..N, so need to index from 1 not 0 */
        ystart[1] = s->water_frac[soil_layer];

        odeint(ystart, N, x1, x2, eps, h1, hmin, &nok, &nbad,
               unsat, drain_layer, p->cond1[soil_layer], p->cond2[soil_layer],
               p->cond3[soil_layer], soil_water_store, rkqs);
        /* ystart is a vector 1..N, so need to index from 1 not 0 */
        new_water_frac = ystart[1];

        /* convert from water fraction to absolute amount (m) */
        change = (s->water_frac[soil_layer] - new_water_frac) * \
                    s->thickness[soil_layer];

        /* update soil layer below with drained liquid */
        f->water_gain[soil_layer+1] += change;
        f->water_loss[soil_layer] += change;
    }

    if (f->water_loss[soil_layer] < 0.0) {
        fprintf(stderr, "waterloss probem in soil_balance\n");
        exit(EXIT_FAILURE);
    }
    free_dvector(ystart, 1, N);

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
    /* determines gravitational water drainage */
    double drainage;

    drainage = calc_soil_conductivity(y[index], cond1, cond2, cond3);

    /* Convert units, soil conductivity is in m s-1 */
    drainage *= SEC_2_HLFHR;

    /* gravitational drainage above field_capacity */
    if (y[index] <= drain_layer) {
        drainage = 0.0;
    }

    /* layer below cannot accept more water than unsat */
    if (drainage > unsat) {
        drainage = unsat;
    }

    /* waterloss from this layer */
    dydt[index] = -drainage;

    return;
}
