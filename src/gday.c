/* ============================================================================
* Generic Decomposition And Yield (GDAY) model.
*
* G'DAY simulates carbon, nutrient and water state and fluxes on either a daily
* sub-daily (30 min) timestep. See below for model
* description.
*
* Paramaeter descriptions are in gday.h
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

#include "gday.h"

int main(int argc, char **argv)
{
    int error = 0;

    /*
     * Setup structures, initialise stuff, e.g. zero fluxes.
     */
    control *c;
    canopy_wk *cw;
    fluxes *f;
    met_arrays *ma;
    met *m;
    params *p;
    state *s;

    c = (control *)malloc(sizeof(control));
    if (c == NULL) {
        fprintf(stderr, "control structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    cw = (canopy_wk *)malloc(sizeof(canopy_wk));
    if (cw == NULL) {
        fprintf(stderr, "canopy wk structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    f = (fluxes *)malloc(sizeof(fluxes));
    if (f == NULL) {
    	fprintf(stderr, "fluxes structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    ma = (met_arrays *)malloc(sizeof(met_arrays));
    if (ma == NULL) {
    	fprintf(stderr, "met arrays structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    m = (met *)malloc(sizeof(met));
    if (m == NULL) {
    	fprintf(stderr, "met structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    p = (params *)malloc(sizeof(params));
    if (p == NULL) {
    	fprintf(stderr, "params structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    s = (state *)malloc(sizeof(state));
    if (s == NULL) {
    	fprintf(stderr, "state structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    initialise_control(c);
    initialise_params(p);
    initialise_fluxes(f);
    initialise_state(s);

    clparser(argc, argv, c);
    /*
     * Read .ini parameter file and meterological data
     */
    error = parse_ini_file(c, p, s);
    if (error != 0) {
        prog_error("Error reading .INI file on line", __LINE__);
    }
    strcpy(c->git_code_ver, build_git_sha);
    if (c->PRINT_GIT) {
        fprintf(stderr, "\n%s\n", c->git_code_ver);
        exit(EXIT_FAILURE);
    }

    /* House keeping! */
    if (c->water_balance == HYDRAULICS && c->sub_daily == FALSE) {
        fprintf(stderr, "You can't run the hydraulics model with daily flag\n");
        exit(EXIT_FAILURE);
    }

    if (c->sub_daily) {
        read_subdaily_met_data(argv, c, ma);
    } else {
        read_daily_met_data(argv, c, ma);
    }

    if (c->spin_up) {
        spin_up_pools(cw, c, f, ma, m, p, s);
    } else {
        run_sim(cw, c, f, ma, m, p, s);
    }

    /* clean up */
    fclose(c->ofp);
    if (c->print_options == SUBDAILY ) {
        fclose(c->ofp_sd);
    }
    fclose(c->ifp);
    if (c->output_ascii == FALSE) {
        fclose(c->ofp_hdr);
    }

    free(cw);
    free(c);
    free(ma->year);
    free(ma->tair);
    free(ma->rain);
    free(ma->tsoil);
    free(ma->co2);
    free(ma->ndep);
    free(ma->wind);
    free(ma->press);
    free(ma->par);
    if (c->sub_daily) {
        free(ma->vpd);
        free(ma->doy);

        /* Clean up hydraulics */
        if (c->water_balance == HYDRAULICS) {
            free(f->soil_conduct);
            free(f->swp);
            free(f->soilR);
            free(f->fraction_uptake);
            free(f->ppt_gain);
            free(f->water_loss);
            free(f->water_gain);

            free(s->water_frac);
            free(s->wetting_bot);
            free(s->wetting_top);
            free(s->thickness);
            free(s->root_mass);
            free(s->root_length);
            free(s->layer_depth);

            free(p->potA);
            free(p->potB);
            free(p->cond1);
            free(p->cond2);
            free(p->cond3);
            free(p->porosity);
            free(p->field_capacity);
        }

    } else {
        free(ma->prjday);
        free(ma->tam);
        free(ma->tpm);
        free(ma->tmin);
        free(ma->tmax);
        free(ma->vpd_am);
        free(ma->vpd_pm);
        free(ma->wind_am);
        free(ma->wind_pm);
        free(ma->par_am);
        free(ma->par_pm);
    }
    free(ma);
    free(m);
    free(p);
    free(s);
    free(f);



    exit(EXIT_SUCCESS);
}

void run_sim(canopy_wk *cw, control *c, fluxes *f, met_arrays *ma, met *m,
             params *p, state *s){

    int    nyr, doy, window_size, i, dummy = 0;
    int    fire_found = FALSE;;
    int    num_disturbance_yrs = 0;

    double fdecay, rdecay, current_limitation, nitfac, year;
    int   *disturbance_yrs = NULL;

    /*
     * potentially allocating 1 extra spot, but will be fine as we always
     * index by num_days
     */
    double *day_length = NULL;
    if ((day_length = (double *)calloc(366, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for day_length\n");
		exit(EXIT_FAILURE);
    }

    if (c->deciduous_model) {
        /* Are we reading in last years average growing season? */
        if (float_eq(s->avg_alleaf, 0.0) &&
            float_eq(s->avg_alstem, 0.0) &&
            float_eq(s->avg_albranch, 0.0) &&
            float_eq(s->avg_alleaf, 0.0) &&
            float_eq(s->avg_alroot, 0.0) &&
            float_eq(s->avg_alcroot, 0.0)) {
            nitfac = 0.0;
            calc_carbon_allocation_fracs(c, f, p, s, nitfac);
        } else {
            f->alleaf = s->avg_alleaf;
            f->alstem = s->avg_alstem;
            f->albranch = s->avg_albranch;
            f->alroot = s->avg_alroot;
            f->alcroot = s->avg_alcroot;
        }
        allocate_stored_c_and_n(f, p, s);
    }

    /* Setup output file */
    if (c->print_options == SUBDAILY && c->spin_up == FALSE) {
        /* open the 30 min outputs file and the daily output files */
        open_output_file(c, c->out_subdaily_fname, &(c->ofp_sd));
        open_output_file(c, c->out_fname, &(c->ofp));

        if (c->output_ascii) {
            write_output_subdaily_header(c, &(c->ofp_sd));
            write_output_header(c, &(c->ofp));
        } else {
            fprintf(stderr, "Nothing implemented for sub-daily binary\n");
            exit(EXIT_FAILURE);
        }
    } else if (c->print_options == DAILY && c->spin_up == FALSE) {
        /* Daily outputs */
        open_output_file(c, c->out_fname, &(c->ofp));

        if (c->output_ascii) {
            write_output_header(c, &(c->ofp));
        } else {
            open_output_file(c, c->out_fname_hdr, &(c->ofp_hdr));
            write_output_header(c, &(c->ofp_hdr));
        }
    } else if (c->print_options == END && c->spin_up == FALSE) {
        /* Final state + param file */
        open_output_file(c, c->out_param_fname, &(c->ofp));
    }

    if (c->water_balance == HYDRAULICS) {
        initialise_roots(f, p, s);
    }

    /*
     * Window size = root lifespan in days...
     * For deciduous species window size is set as the length of the
     * growing season in the main part of the code
     */
    window_size = (int)(1.0 / p->rdecay * NDAYS_IN_YR);
    sma_obj *hw = sma(SMA_NEW, window_size).handle;
    if (s->prev_sma > -900) {
        for (i = 0; i < window_size; i++) {
            sma(SMA_ADD, hw, s->prev_sma);
        }
    }
    /* Set up SMA
     *  - If we don't have any information about the N & water limitation, i.e.
     *    as would be the case with spin-up, assume that there is no limitation
     *    to begin with.
     */
    if (s->prev_sma < -900)
        s->prev_sma = 1.0;

    /*
     * Params are defined in per year, needs to be per day. Important this is
     * done here as rate constants elsewhere in the code are assumed to be in
     * units of days not years
     */
    correct_rate_constants(p, FALSE);
    day_end_calculations(c, p, s, -99, TRUE);

    if (c->sub_daily) {
        initialise_soils_sub_daily(c, f, p, s);
    } else {
        initialise_soils_day(c, f, p, s);
    }

    if (c->water_balance == HYDRAULICS) {
        double root_zone_total, water_content;

        // Update the soil water storage
        root_zone_total = 0.0;
        for (i = 0; i < p->n_layers; i++) {

            // water content of soil layer (m)
            water_content = s->water_frac[i] * s->thickness[i];

            // update old GDAY effective two-layer buckets
            // - this is just for outputting, these aren't used.
            if (i == 0) {
                s->pawater_topsoil = water_content * M_TO_MM;
            } else {
                root_zone_total += water_content * M_TO_MM;
            }
        }
        s->pawater_root = root_zone_total;
    } else {
        s->pawater_root = p->wcapac_root;
        s->pawater_topsoil = p->wcapac_topsoil;
    }

    if (c->fixed_lai) {
        s->lai = p->fix_lai;
    } else {
        s->lai = MAX(0.01, (p->sla * M2_AS_HA / KG_AS_TONNES /
                            p->cfracts * s->shoot));
    }

    if (c->disturbance) {
        if ((disturbance_yrs = (int *)calloc(1, sizeof(int))) == NULL) {
            fprintf(stderr,"Error allocating space for disturbance_yrs\n");
    		exit(EXIT_FAILURE);
        }
        figure_out_years_with_disturbances(c, ma, p, &disturbance_yrs,
                                           &num_disturbance_yrs);
    }

    /* ====================== **
    **   Y E A R    L O O P   **
    ** ====================== */
    c->day_idx = 0;
    c->hour_idx = 0;



    for (nyr = 0; nyr < c->num_years; nyr++) {

        if (c->sub_daily) {
            year = ma->year[c->hour_idx];
        } else {
            year = ma->year[c->day_idx];
        }
        if (is_leap_year(year))
            c->num_days = 366;
        else
            c->num_days = 365;

        calculate_daylength(c->num_days, p->latitude, *(&day_length));

        if (c->deciduous_model) {
            phenology(c, f, ma, p, s, day_length);

            /* Change window size to length of growing season */
            sma(SMA_FREE, hw);
            hw = sma(SMA_NEW, p->growing_seas_len).handle;
            if (s->prev_sma > -900) {
                for (i = 0; i < p->growing_seas_len; i++) {
                    sma(SMA_ADD, hw, s->prev_sma);
                }
            }

            zero_stuff(c, s);
        }
        /* =================== **
        **   D A Y   L O O P   **
        ** =================== */
        for (doy = 0; doy < c->num_days; doy++) {

            //if (year == 2002 && doy+1 == 240) {
            //    c->pdebug = TRUE;
            //}
            if (! c->sub_daily) {
                unpack_met_data(c, f, ma, m, dummy, day_length[doy]);
            }
            calculate_litterfall(c, f, p, s, doy, &fdecay, &rdecay);

            if (c->disturbance && p->disturbance_doy == doy+1) {
                /* Fire Disturbance? */
                fire_found = FALSE;
                fire_found = check_for_fire(c, f, p, s, year, disturbance_yrs,
                                            num_disturbance_yrs);

                if (fire_found) {
                    fire(c, f, p, s);
                    /*
                     * This will only work for evergreen, but that is fine
                     * this should be removed after KSCO is done
                     */
                    sma(SMA_FREE, hw);
                    hw = sma(SMA_NEW, window_size).handle;
                    if (s->prev_sma > -900) {
                        for (i = 0; i < window_size; i++) {
                            sma(SMA_ADD, hw, s->prev_sma);
                        }
                    }
                }
            } else if (c->hurricane &&
                p->hurricane_yr == year &&
                p->hurricane_doy == doy) {

                /* Hurricane? */
                hurricane(f, p, s);
            }


            calc_day_growth(cw, c, f, ma, m, p, s, day_length[doy],
                            doy, fdecay, rdecay);

            /*printf("%d %f %f\n", doy, f->gpp*100, s->lai);*/
            calculate_csoil_flows(c, f, p, s, m->tsoil, doy);
            calculate_nsoil_flows(c, f, p, s, doy);

            /* update stress SMA */
            if (c->deciduous_model && s->leaf_out_days[doy] > 0.0) {
                 /*
                  * Allocation is annually for deciduous "tree" model, but we
                  * need to keep a check on stresses during the growing season
                  * and the LAI figure out limitations during leaf growth period.
                  * This also applies for deciduous grasses, need to do the
                  * growth stress calc for grasses here too.
                  */
                current_limitation = calculate_growth_stress_limitation(p, s);
                sma(SMA_ADD, hw, current_limitation);
                s->prev_sma = sma(SMA_MEAN, hw).sma;
            } else if (c->deciduous_model == FALSE) {
                current_limitation = calculate_growth_stress_limitation(p, s);
                sma(SMA_ADD, hw, current_limitation);
                s->prev_sma = sma(SMA_MEAN, hw).sma;
            }

            /*
             * if grazing took place need to reset "stress" running mean
             * calculation for grasses
             */
            if (c->grazing == 2 && p->disturbance_doy == doy+1) {
                sma(SMA_FREE, hw);
                hw = sma(SMA_NEW, p->growing_seas_len).handle;
            }

            /* Turn off all N calculations */
            if (c->ncycle == FALSE)
                reset_all_n_pools_and_fluxes(f, s);

            /* calculate C:N ratios and increment annual flux sum */
            day_end_calculations(c, p, s, c->num_days, FALSE);

            if (c->print_options == SUBDAILY && c->spin_up == FALSE) {
                write_daily_outputs_ascii(c, f, s, year, doy+1);
            } else if (c->print_options == DAILY && c->spin_up == FALSE) {
                if(c->output_ascii)
                    write_daily_outputs_ascii(c, f, s, year, doy+1);
                else
                    write_daily_outputs_binary(c, f, s, year, doy+1);
            }
            c->day_idx++;


            printf("%d %d %f", (int)year, doy, s->water_frac[0] * s->thickness[0] * M_TO_MM);
            //printf("%d %d %f", (int)year, doy, s->water_frac[0]);
            for (i = 1; i < p->n_layers; i++) {
            //
            //    //printf(" %f", s->water_frac[i] * s->thickness[i] * M_TO_MM);
                printf(" %f", s->water_frac[i]);
            //
            }
            //printf("\n");
            //printf("%d %d %lf %lf %lf\n", (int)year, doy, s->saved_swp, s->wtfac_root, f->gpp*100);

            //printf("%d %d %lf %lf %lf %lf\n", (int)year, doy, f->gpp*100, f->transpiration, s->wtfac_root, s->saved_swp);
            //printf("%d %d %lf %lf %lf\n", (int)year, doy, f->gpp*100, f->transpiration, s->wtfac_root);

            printf("\n");
            /* ======================= **
            **   E N D   O F   D A Y   **
            ** ======================= */
        }

        /* Allocate stored C&N for the following year */
        if (c->deciduous_model) {
            calculate_average_alloc_fractions(f, s, p->growing_seas_len);
            allocate_stored_c_and_n(f, p, s);
        }

        // Adjust rooting distribution at the end of the year to account for
        // growth of new roots. It is debatable when this should be done. I've
        // picked the year end for computation reasons and probably because
        // plants wouldn't do this as dynamcially as on a daily basis. Probably
        if (c->water_balance == HYDRAULICS) {
            update_roots(c, p, s);
        }
    }

    /* ========================= **
    **   E N D   O F   Y E A R   **
    ** ========================= */
    correct_rate_constants(p, TRUE);

    if (c->print_options == END && c->spin_up == FALSE) {
        write_final_state(c, p, s);
    }

    sma(SMA_FREE, hw);
    free(day_length);
    if (c->disturbance) {
        free(disturbance_yrs);
    }

    return;


}

void spin_up_pools(canopy_wk *cw, control *c, fluxes *f, met_arrays *ma, met *m,
                   params *p, state *s){
    /* Spin up model plant & soil pools to equilibrium.

    - Examine sequences of 50 years and check if C pools are changing
      by more than 0.005 units per 1000 yrs. Note this check is done in
      units of: kg m-2.

    References:
    ----------
    Adapted from...
    * Murty, D and McMurtrie, R. E. (2000) Ecological Modelling, 134,
      185-205, specifically page 196.
    */
    double tol = 5E-03;
    double prev_plantc = 99999.9;
    double prev_soilc = 99999.9;
    int i, cntrl_flag;
    /* check for convergences in units of kg/m2 */
    double conv = TONNES_HA_2_KG_M2;

    /* Final state + param file */
    open_output_file(c, c->out_param_fname, &(c->ofp));

    /* If we are prescribing disturbance, first allow the forest to establish */
    if (c->disturbance) {
        cntrl_flag = c->disturbance;
        c->disturbance = FALSE;
        /*  200 years (50 yrs x 4 cycles) */
        for (i = 0; i < 4; i++) {
            run_sim(cw, c, f, ma, m, p, s); /* run GDAY */
        }
        c->disturbance = cntrl_flag;
    }

    fprintf(stderr, "Spinning up the model...\n");
    while (TRUE) {
        if (fabs((prev_plantc*conv) - (s->plantc*conv)) < tol &&
            fabs((prev_soilc*conv) - (s->soilc*conv)) < tol) {
            break;
        } else {
            prev_plantc = s->plantc;
            prev_soilc = s->soilc;

            /* 1000 years (50 yrs x 20 cycles) */
            for (i = 0; i < 20; i++) {
                run_sim(cw, c, f, ma, m, p, s); /* run GDAY */
            }

            /* Have we reached a steady state? */
            fprintf(stderr,
              "Spinup: Plant C - %f, Soil C - %f\n", s->plantc, s->soilc);
        }
    }
    write_final_state(c, p, s);

    return;
}

void clparser(int argc, char **argv, control *c) {
    int i;

    for (i = 1; i < argc; i++) {
        if (*argv[i] == '-') {
            if (!strncasecmp(argv[i], "-p", 2)) {
			    strcpy(c->cfg_fname, argv[++i]);
            } else if (!strncasecmp(argv[i], "-s", 2)) {
                c->spin_up = TRUE;
            } else if (!strncasecmp(argv[i], "-ver", 4)) {
                c->PRINT_GIT = TRUE;
            } else if (!strncasecmp(argv[i], "-u", 2) ||
                       !strncasecmp(argv[i], "-h", 2)) {
                usage(argv);
                exit(EXIT_FAILURE);
            } else {
                fprintf(stderr, "%s: unknown argument on command line: %s\n",
                               argv[0], argv[i]);
                usage(argv);
                exit(EXIT_FAILURE);
            }
        }
    }
    return;
}


void usage(char **argv) {
    fprintf(stderr, "\n========\n");
    fprintf(stderr, " USAGE:\n");
    fprintf(stderr, "========\n");
    fprintf(stderr, "%s [options]\n", argv[0]);
    fprintf(stderr, "\n\nExpected input file is a .ini/.cfg style param file, passed with the -p flag .\n");
    fprintf(stderr, "\nThe options are:\n");
    fprintf(stderr, "\n++General options:\n" );
    fprintf(stderr, "[-ver          \t] Print the git hash tag.]\n");
    fprintf(stderr, "[-p       fname\t] Location of parameter file (.ini/.cfg).]\n");
    fprintf(stderr, "[-s            \t] Spin-up GDAY, when it the model is finished it will print the final state to the param file.]\n");
    fprintf(stderr, "\n++Print this message:\n" );
    fprintf(stderr, "[-u/-h         \t] usage/help]\n");

    return;
}





void correct_rate_constants(params *p, int output) {
    /* adjust rate constants for the number of days in years */

    if (output) {
        p->rateuptake *= NDAYS_IN_YR;
        p->rateloss *= NDAYS_IN_YR;
        p->retransmob *= NDAYS_IN_YR;
        p->fdecay *= NDAYS_IN_YR;
        p->fdecaydry *= NDAYS_IN_YR;
        p->crdecay *= NDAYS_IN_YR;
        p->rdecay *= NDAYS_IN_YR;
        p->rdecaydry *= NDAYS_IN_YR;
        p->bdecay *= NDAYS_IN_YR;
        p->wdecay *= NDAYS_IN_YR;
        p->sapturnover *= NDAYS_IN_YR;
        p->kdec1 *= NDAYS_IN_YR;
        p->kdec2 *= NDAYS_IN_YR;
        p->kdec3 *= NDAYS_IN_YR;
        p->kdec4 *= NDAYS_IN_YR;
        p->kdec5 *= NDAYS_IN_YR;
        p->kdec6 *= NDAYS_IN_YR;
        p->kdec7 *= NDAYS_IN_YR;
        p->nuptakez *= NDAYS_IN_YR;
        p->nmax *= NDAYS_IN_YR;
    } else {
        p->rateuptake /= NDAYS_IN_YR;
        p->rateloss /= NDAYS_IN_YR;
        p->retransmob /= NDAYS_IN_YR;
        p->fdecay /= NDAYS_IN_YR;
        p->fdecaydry /= NDAYS_IN_YR;
        p->crdecay /= NDAYS_IN_YR;
        p->rdecay /= NDAYS_IN_YR;
        p->rdecaydry /= NDAYS_IN_YR;
        p->bdecay /= NDAYS_IN_YR;
        p->wdecay /= NDAYS_IN_YR;
        p->sapturnover /= NDAYS_IN_YR;
        p->kdec1 /= NDAYS_IN_YR;
        p->kdec2 /= NDAYS_IN_YR;
        p->kdec3 /= NDAYS_IN_YR;
        p->kdec4 /= NDAYS_IN_YR;
        p->kdec5 /= NDAYS_IN_YR;
        p->kdec6 /= NDAYS_IN_YR;
        p->kdec7 /= NDAYS_IN_YR;
        p->nuptakez /= NDAYS_IN_YR;
        p->nmax /= NDAYS_IN_YR;
    }

    return;
}


void reset_all_n_pools_and_fluxes(fluxes *f, state *s) {
    /*
        If the N-Cycle is turned off the way I am implementing this is to
        do all the calculations and then reset everything at the end. This is
        a waste of resources but saves on multiple IF statements.
    */

    /*
    ** State
    */
    s->shootn = 0.0;
    s->rootn = 0.0;
    s->crootn = 0.0;
    s->branchn = 0.0;
    s->stemnimm = 0.0;
    s->stemnmob = 0.0;
    s->structsurfn = 0.0;
    s->metabsurfn = 0.0;
    s->structsoiln = 0.0;
    s->metabsoiln = 0.0;
    s->activesoiln = 0.0;
    s->slowsoiln = 0.0;
    s->passivesoiln = 0.0;
    s->inorgn = 0.0;
    s->stemn = 0.0;
    s->stemnimm = 0.0;
    s->stemnmob = 0.0;
    s->nstore = 0.0;

    /*
    ** Fluxes
    */
    f->nuptake = 0.0;
    f->nloss = 0.0;
    f->npassive = 0.0;
    f->ngross = 0.0;
    f->nimmob = 0.0;
    f->nlittrelease = 0.0;
    f->nmineralisation = 0.0;
    f->npleaf = 0.0;
    f->nproot = 0.0;
    f->npcroot = 0.0;
    f->npbranch = 0.0;
    f->npstemimm = 0.0;
    f->npstemmob = 0.0;
    f->deadleafn = 0.0;
    f->deadrootn = 0.0;
    f->deadcrootn = 0.0;
    f->deadbranchn = 0.0;
    f->deadstemn = 0.0;
    f->neaten = 0.0;
    f->nurine = 0.0;
    f->leafretransn = 0.0;
    f->n_surf_struct_litter = 0.0;
    f->n_surf_metab_litter = 0.0;
    f->n_soil_struct_litter = 0.0;
    f->n_soil_metab_litter = 0.0;
    f->n_surf_struct_to_slow = 0.0;
    f->n_soil_struct_to_slow = 0.0;
    f->n_surf_struct_to_active = 0.0;
    f->n_soil_struct_to_active = 0.0;
    f->n_surf_metab_to_active = 0.0;
    f->n_surf_metab_to_active = 0.0;
    f->n_active_to_slow = 0.0;
    f->n_active_to_passive = 0.0;
    f->n_slow_to_active = 0.0;
    f->n_slow_to_passive = 0.0;
    f->n_passive_to_active = 0.0;

    return;
}

void zero_stuff(control *c, state *s) {
    s->shoot = 0.0;
    s->shootn = 0.0;
    s->shootnc = 0.0;
    s->lai = 0.0;
    s->cstore = 0.0;
    s->nstore = 0.0;
    s->anpp = 0.0;

    if (c->deciduous_model) {
        s->avg_alleaf = 0.0;
        s->avg_alroot = 0.0;
        s->avg_alcroot = 0.0;
        s->avg_albranch  = 0.0;
        s->avg_alstem = 0.0;
    }
    return;
}

void day_end_calculations(control *c, params *p, state *s, int days_in_year,
                          int init) {
    /* Calculate derived values from state variables.

    Parameters:
    -----------
    day : integer
        day of simulation

    INIT : logical
        logical defining whether it is the first day of the simulation
    */

    /* update N:C of plant pool */
    if (float_eq(s->shoot, 0.0))
        s->shootnc = 0.0;
    else
        s->shootnc = s->shootn / s->shoot;

    /* Explicitly set the shoot N:C */
    if (c->ncycle == FALSE)
        s->shootnc = p->prescribed_leaf_NC;

    if (float_eq(s->root, 0.0))
        s->rootnc = 0.0;
    else
        s->rootnc = MAX(0.0, s->rootn / s->root);

    /* total plant, soil & litter nitrogen */
    s->soiln = s->inorgn + s->activesoiln + s->slowsoiln + s->passivesoiln;
    s->litternag = s->structsurfn + s->metabsurfn;
    s->litternbg = s->structsoiln + s->metabsoiln;
    s->littern = s->litternag + s->litternbg;
    s->plantn = s->shootn + s->rootn + s->crootn + s->branchn + s->stemn;
    s->totaln = s->plantn + s->littern + s->soiln;

    /* total plant, soil, litter and system carbon */
    s->soilc = s->activesoil + s->slowsoil + s->passivesoil;
    s->littercag = s->structsurf + s->metabsurf;
    s->littercbg = s->structsoil + s->metabsoil;
    s->litterc = s->littercag + s->littercbg;
    s->plantc = s->root + s->croot + s->shoot + s->stem + s->branch;
    s->totalc = s->soilc + s->litterc + s->plantc;

    /* optional constant passive pool */
    if (c->passiveconst) {
        s->passivesoil = p->passivesoilz;
        s->passivesoiln = p->passivesoilnz;
    }

    if (init == FALSE)
        /* Required so max leaf & root N:C can depend on Age */
        s->age += 1.0 / days_in_year;

    return;
}

void unpack_met_data(control *c, fluxes *f, met_arrays *ma, met *m, int hod,
                     double day_length) {

    double c1, c2;

    /* unpack met forcing */
    if (c->sub_daily) {
        m->rain = ma->rain[c->hour_idx];
        m->wind = ma->wind[c->hour_idx];
        m->press = ma->press[c->hour_idx] * KPA_2_PA;
        m->vpd = ma->vpd[c->hour_idx] * KPA_2_PA;
        m->tair = ma->tair[c->hour_idx];
        m->tsoil = ma->tsoil[c->hour_idx];
        m->par = ma->par[c->hour_idx];
        m->sw_rad = ma->par[c->hour_idx] * PAR_2_SW; /* W m-2 */
        m->Ca = ma->co2[c->hour_idx];

        /* NDEP is per 30 min so need to sum 30 min data */
        if (hod == 0) {
            m->ndep = ma->ndep[c->hour_idx];
            m->nfix = ma->nfix[c->hour_idx];
        } else {
            m->ndep += ma->ndep[c->hour_idx];
            m->nfix += ma->nfix[c->hour_idx];
        }
    } else {
        m->Ca = ma->co2[c->day_idx];
        m->tair = ma->tair[c->day_idx];
        m->tair_am = ma->tam[c->day_idx];
        m->tair_pm = ma->tpm[c->day_idx];
        m->par = ma->par_am[c->day_idx] + ma->par_pm[c->day_idx];

        /* Conversion factor for PAR to SW rad */
        c1 = MJ_TO_J * J_2_UMOL / (day_length * 60.0 * 60.0) * PAR_2_SW;
        c2 = MJ_TO_J * J_2_UMOL / (day_length / 2.0 * 60.0 * 60.0) * PAR_2_SW;
        m->sw_rad = m->par * c1;
        m->sw_rad_am = ma->par_am[c->day_idx] * c2;
        m->sw_rad_pm = ma->par_pm[c->day_idx] * c2;
        m->rain = ma->rain[c->day_idx];
        m->vpd_am = ma->vpd_am[c->day_idx] * KPA_2_PA;
        m->vpd_pm = ma->vpd_pm[c->day_idx] * KPA_2_PA;
        m->wind_am = ma->wind_am[c->day_idx];
        m->wind_pm = ma->wind_pm[c->day_idx];
        m->press = ma->press[c->day_idx] * KPA_2_PA;
        m->ndep = ma->ndep[c->day_idx];
        m->nfix = ma->nfix[c->day_idx];
        m->tsoil = ma->tsoil[c->day_idx];
        m->Tk_am = ma->tam[c->day_idx] + DEG_TO_KELVIN;
        m->Tk_pm = ma->tpm[c->day_idx] + DEG_TO_KELVIN;

        /*printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
               m->Ca, m->tair, m->tair_am, m->tair_pm, m->par, m->sw_rad,
               m->sw_rad_am, m->sw_rad_pm, m->rain, m->vpd_am, m->vpd_pm,
               m->wind_am, m->wind_pm, m->press, m->ndep, m->tsoil, m->Tk_am,
               m->Tk_pm);*/

    }

    /* N deposition + biological N fixation */
    f->ninflow = m->ndep + m->nfix;

    return;
}
