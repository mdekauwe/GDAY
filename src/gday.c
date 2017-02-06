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
    nrutil *nr;
    fast_spinup *fs;

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

    nr = (nrutil *)malloc(sizeof(nrutil));
    if (nr == NULL) {
        fprintf(stderr, "nrutil structure: Not allocated enough memory!\n");
        exit(EXIT_FAILURE);
    }

    fs = (fast_spinup *)malloc(sizeof(fast_spinup));
    if (cw == NULL) {
        fprintf(stderr, "fast spinup structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    // potentially allocating 1 extra spot, but will be fine as we always
    // index by num_days
    if ((s->day_length = (double *)calloc(366, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for day_length\n");
		exit(EXIT_FAILURE);
    }

    initialise_control(c);
    initialise_params(p);
    initialise_fluxes(f);
    initialise_state(s);
    initialise_nrutil(nr);

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

    if (c->water_balance == HYDRAULICS) {
        allocate_numerical_libs_stuff(nr);
        initialise_roots(f, p, s);
        setup_hydraulics_arrays(f, p, s);
    }

    if (c->sub_daily) {
        read_subdaily_met_data(argv, c, ma);
        fill_up_solar_arrays(cw, c, ma, p);
    } else {
        read_daily_met_data(argv, c, ma);
    }


    if (c->spin_up) {
        spin_up_pools(cw, c, f, fs, ma, m, p, s, nr);
    } else {
        run_sim(cw, c, f, fs, ma, m, p, s, nr);
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
        free(cw->cz_store);
        free(cw->ele_store);
        free(cw->df_store);

        /* Clean up hydraulics */
        if (c->water_balance == HYDRAULICS) {
            free(f->soil_conduct);
            free(f->swp);
            free(f->soilR);
            free(f->fraction_uptake);
            free(f->ppt_gain);
            free(f->water_loss);
            free(f->water_gain);
            free(f->est_evap);
            free(s->water_frac);
            free(s->wetting_bot);
            free(s->wetting_top);
            free(p->potA);
            free(p->potB);
            free(p->cond1);
            free(p->cond2);
            free(p->cond3);
            free(p->porosity);
            free(p->field_capacity);
            free(s->thickness);
            free(s->root_mass);
            free(s->root_length);
            free(s->layer_depth);


            free_dvector(nr->y, 1, nr->N);
            free_dvector(nr->ystart, 1, nr->N);
            free_dvector(nr->dydx, 1, nr->N);
			free_dvector(nr->yscal, 1, nr->N);
            free_dvector(nr->xp, 1, nr->kmax);
            free_dmatrix(nr->yp, 1, nr->N, 1, nr->kmax);
            free_dvector(nr->ytemp, 1, nr->N);
        	free_dvector(nr->ak6, 1, nr->N);
        	free_dvector(nr->ak5, 1, nr->N);
        	free_dvector(nr->ak4, 1, nr->N);
        	free_dvector(nr->ak3, 1, nr->N);
        	free_dvector(nr->ak2, 1, nr->N);
            free_dvector(nr->yerr, 1, nr->N);
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
    free(s->day_length);
    free(ma);
    free(m);
    free(p);
    free(s);
    free(f);
    free(fs);

    exit(EXIT_SUCCESS);
}

void run_sim(canopy_wk *cw, control *c, fluxes *f, fast_spinup *fs,
             met_arrays *ma, met *m, params *p, state *s, nrutil *nr) {

    int    nyr, doy, window_size, i, dummy = 0;
    int    fire_found = FALSE;;
    int    num_disturbance_yrs = 0;

    double fdecay, rdecay, current_limitation, nitfac, year;
    int   *disturbance_yrs = NULL;

    if (c->deciduous_model) {
        /* Are we reading in last years average growing season? */
        if (float_eq(s->avg_alleaf, 0.0) &&
            float_eq(s->avg_alstem, 0.0) &&
            float_eq(s->avg_albranch, 0.0) &&
            float_eq(s->avg_alleaf, 0.0) &&
            float_eq(s->avg_alroot, 0.0) &&
            float_eq(s->avg_alcroot, 0.0)) {
            nitfac = 0.0;
            calc_carbon_allocation_fracs(c, f, fs, p, s, nitfac);
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

        calculate_daylength(s, c->num_days, p->latitude);

        if (c->deciduous_model) {
            phenology(c, f, ma, p, s);

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

            //if (year == 2001 && doy+1 == 230) {
            //    c->pdebug = TRUE;
            //}


            if (! c->sub_daily) {
                unpack_met_data(c, f, ma, m, dummy, s->day_length[doy]);
            }
            calculate_litterfall(c, f, fs, p, s, doy, &fdecay, &rdecay);

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


            calc_day_growth(cw, c, f, fs, ma, m, nr, p, s, s->day_length[doy],
                            doy, fdecay, rdecay);

            //printf("%d %f %f\n", doy, f->gpp*100, s->lai);
            calculate_csoil_flows(c, f, fs, p, s, m->tsoil, doy);
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

            // Step 2: Store the time-varying variables
            if (c->spinup_method == SAS) {
                fs->npp_ss += f->npp;
                fs->ndays ++;
                fs->shoot_nc += s->shootn / s->shoot;
                fs->root_nc += s->rootn / s->root;
                fs->branch_nc += s->branchn / s->branch;
                if (s->croot > 0.0) {
                    fs->croot_nc += s->crootn / s->croot;
                } else {
                    fs->croot_nc = 0.0;
                }
                fs->stem_nc += s->stemn / s->stem;
                if (s->stemnmob > 0.0) {
                    fs->stemnmob_ratio += s->stemnmob / s->stem;
                } else {
                    fs->stemnmob_ratio = 0.0;
                }
                if (s->stemnimm > 0.0) {
                    fs->stemnimm_ratio += s->stemnimm / s->stem;
                } else {
                    fs->stemnimm_ratio = 0.0;
                }

                if (s->metabsoil > 0.0) {
                    fs->metablsoil_nc += s->metabsoiln / s->metabsoil;
                } else {
                    fs->metablsoil_nc += 0.0;
                }

                if (s->metabsurf > 0.0) {
                    fs->metabsurf_nc += s->metabsurfn / s->metabsurf;
                } else {
                    fs->metabsurf_nc += 0.0;
                }

                fs->structsoil_nc += s->structsoiln / s->structsoil;
                fs->structsurf_nc += s->structsurfn / s->structsurf;
                fs->activesoil_nc += s->activesoiln / s->activesoil;
                fs->slowsoil_nc += s->slowsoiln / s->slowsoil;
                fs->passivesoil_nc += s->passivesoiln / s->passivesoil;
            }






            c->day_idx++;


            //printf("%d %d %f", (int)year, doy, s->water_frac[0] * s->thickness[0] * M_TO_MM);
            //printf("%d %d %f", (int)year, doy, s->water_frac[0]);
            //for (i = 1; i < p->n_layers; i++) {
            //
            //    //printf(" %f", s->water_frac[i] * s->thickness[i] * M_TO_MM);
            //    printf(" %f", s->water_frac[i]);
            //
            //}
            //printf("\n");
            //printf("%d %d %lf %lf %lf\n", (int)year, doy, s->saved_swp, s->wtfac_root, f->gpp*100);

            //printf("%d %d %lf %lf %lf %lf\n", (int)year, doy, f->gpp*100, f->transpiration, s->wtfac_root, s->saved_swp);
            //printf("%d %d %lf %lf %lf\n", (int)year, doy, f->gpp*100, f->transpiration, s->wtfac_root);


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
    if (c->disturbance) {
        free(disturbance_yrs);
    }

    return;


}

void spin_up_pools(canopy_wk *cw, control *c, fluxes *f, fast_spinup *fs,
                   met_arrays *ma, met *m, params *p, state *s, nrutil *nr) {
    /* Spin up model plant & soil pools to equilibrium.

    - Examine sequences of 50 years and check if C pools are changing
      by more than 0.005 units per 1000 yrs.

    References:
    ----------
    Adapted from...
    * Murty, D and McMurtrie, R. E. (2000) Ecological Modelling, 134,
      185-205, specifically page 196.
    */
    double tol = 5E-03;
    double prev_plantc = 99999.9;
    double prev_soilc = 99999.9;
    int    i, cntrl_flag;

    /* Final state + param file */
    open_output_file(c, c->out_param_fname, &(c->ofp));

    /* If we are prescribing disturbance, first allow the forest to establish */
    if (c->disturbance) {
        cntrl_flag = c->disturbance;
        c->disturbance = FALSE;
        /*  200 years (50 yrs x 4 cycles) */
        for (i = 0; i < 4; i++) {
            run_sim(cw, c, f, fs, ma, m, p, s, nr); /* run GDAY */
        }
        c->disturbance = cntrl_flag;
    }

    fprintf(stderr, "Spinning up the model...\n");

    if (c->spinup_method == BRUTE) {
        printf("HERE %d\n", c->spinup_method);
    } else {
        printf("NO HERE\n");
    }

    exit(1);


    if (c->spinup_method == BRUTE) {

        while (TRUE) {
            if (fabs(prev_plantc - s->plantc) < tol &&
                fabs(prev_soilc - s->soilc) < tol) {
                break;
            } else {
                prev_plantc = s->plantc;
                prev_soilc = s->soilc;

                /* 1000 years (50 yrs x 20 cycles) */
                for (i = 0; i < 20; i++) {
                    run_sim(cw, c, f, fs, ma, m, p, s, nr); /* run GDAY */
                }

                /* Have we reached a steady state? */
                fprintf(stderr,
                  "Spinup: Plant C - %f, Soil C - %f\n", s->plantc, s->soilc);
            }

            /* total plant, soil, litter and system carbon */
            s->soilc = s->activesoil + s->slowsoil + s->passivesoil;
            s->littercag = s->structsurf + s->metabsurf;
            s->littercbg = s->structsoil + s->metabsoil;
            s->litterc = s->littercag + s->littercbg;
            s->plantc = s->root + s->croot + s->shoot + s->stem + s->branch;
            s->totalc = s->soilc + s->litterc + s->plantc;

        }

    } else if (c->spinup_method == SAS) {
        //
        // Semi-analytical solution (SAS) to accelerate model spin-up of
        // carbon–nitrogen pools, following Xia et al. (2013) GMD.
        //
        sas_spinup(cw, c, f, fs, ma, m, p, s, nr);
    }
    exit(1);

    write_final_state(c, p, s);

    return;
}

void sas_spinup(canopy_wk *cw, control *c, fluxes *f, fast_spinup *fs,
                met_arrays *ma, met *m, params *p, state *s, nrutil *nr) {
    //
    // Semi-analytical solution (SAS) to accelerate model spin-up of
    // carbon–nitrogen pools, following Xia et al. (2013) GMD.
    //

    double cleaf0, cwood0, croot0, criteria, arg1, arg2, arg3;
    double NPP, mu_af, mu_ar, mu_acr, mu_ab, mu_aw, mu_lf, mu_lr, mu_lcr;
    double mu_lb, mu_lw, shootX, rootX, crootX, branchX, stemX, wood, woodX;
    double mu_ass1, mu_ass2, mu_ass3, leaf_material, wood_material, mu_as1;
    double surf_struct_litter, structout_surf, structout_soil;
    double surf_struct_to_slow, surf_struct_to_active;
    double soil_struct_to_slow, soil_struct_litter;
    double soil_metab_litter, metabsurfX, metabsoilX;
    double structsurfX, structsoilX, surf_metab_to_active, soil_metab_to_active;
    double co2_to_air0, co2_to_air1, co2_to_air2, co2_to_air3, co2_to_air4;
    double co2_to_air5, co2_to_air6;
    double activesoilX, slowsoilX, passivesoilX, passive_to_active;
    double c_into_active, slow_to_active, slow_to_passive, slowout;
    double activeout, frac_microb_resp, c_into_passive;
    double active_to_slow, active_to_passive, c_into_slow, mu_fmleaf, mu_fmroot;
    double leafgrowth, rootgrowth, crootgrowth, branchgrowth, stemgrowth;
    double deadleaves, deadroots, deadcroots, deadbranches, deadstems;
    double mu_decayrate0, mu_decayrate1, mu_decayrate2, mu_decayrate3, mu_decayrate4;
    double mu_decayrate5, mu_decayrate6, surf_metab_litter, soil_struct_to_active;
    double total_days, deadsapwood, sapwoodX, new_passive;
    double prev_passivec = 99999.9;
    int    i, cntrl_flag;

    // Step 1: Initial spin
    // - we first need to achieve steady state plant pools (or NPP is an
    //   alternative, I didn't test that).
    zero_fast_spinup_stuff(fs);
    run_sim(cw, c, f, fs, ma, m, p, s, nr); /* run GDAY */
    cleaf0 = s->shoot;
    cwood0 = s->branch + s->croot + s->stem;
    croot0 = s->root;
    criteria = 0.01 * (cleaf0 + cwood0 + croot0);

    while (TRUE) {
        zero_fast_spinup_stuff(fs);
        run_sim(cw, c, f, fs, ma, m, p, s, nr); /* run GDAY */
        arg1 = fabs((s->shoot - cleaf0) / s->shoot);
        wood = s->branch + s->croot + s->stem;
        arg2 = fabs((wood - cwood0) / wood);
        arg3 = fabs((s->root - croot0) / s->root);

        if ((arg1 + arg2 + arg3) < criteria) {
            break;
        } else {
            cleaf0 = s->shoot;
            cwood0 = s->branch + s->croot + s->stem;
            croot0 = s->root;
            criteria = 0.01 * (cleaf0 + cwood0 + croot0);
        }
    }

    // Step 2: store the time varying vars is done internal in the rest of the
    //         code.

    // Step 3: Analytically solve the steady state pools

    // Calculate the mean time-varying variables
    total_days = (double)fs->ndays;
    mu_af = fs->alloc[AF] / total_days;
    mu_lf = fs->loss[LF] / total_days;
    mu_ar = fs->alloc[AR] / total_days;
    mu_lr = fs->loss[LR] / total_days;
    mu_acr = fs->alloc[ACR] / total_days;
    mu_lcr = fs->loss[LCR] / total_days;
    mu_ab = fs->alloc[AB] / total_days;
    mu_lb = fs->loss[LB] / total_days;
    mu_aw = fs->alloc[AW] / total_days;
    mu_lw = fs->loss[LW] / total_days;
    mu_decayrate0 = fs->dr[0] / total_days;
    mu_decayrate1 = fs->dr[1] / total_days;
    mu_decayrate2 = fs->dr[2] / total_days;
    mu_decayrate3 = fs->dr[3] / total_days;
    mu_decayrate4 = fs->dr[4] / total_days;
    mu_decayrate5 = fs->dr[5] / total_days;
    mu_decayrate6 = fs->dr[5] / total_days;
    mu_fmleaf = fs->alloc[S1] / total_days;
    mu_fmroot = fs->alloc[S2] / total_days;

    fs->shoot_nc /= total_days;
    fs->root_nc /= total_days;
    fs->branch_nc /= total_days;
    fs->croot_nc /= total_days;
    fs->stem_nc /= total_days;
    fs->stemnmob_ratio /= total_days;
    fs->stemnimm_ratio /= total_days;
    fs->metablsoil_nc /= total_days;
    fs->metabsurf_nc /= total_days;
    fs->structsoil_nc /= total_days;
    fs->structsurf_nc /= total_days;
    fs->activesoil_nc /= total_days;
    fs->slowsoil_nc /= total_days;
    fs->passivesoil_nc /= total_days;

    // Steady-state NPP
    NPP = fs->npp_ss / total_days;

    // Solve the C pools
    leafgrowth = (NPP * mu_af);
    deadleaves = (s->shoot * mu_lf);
    shootX = leafgrowth - deadleaves;

    rootgrowth = (NPP * mu_ar);
    deadroots = (s->root * mu_lr);
    rootX = rootgrowth - deadroots;

    crootgrowth = (NPP * mu_acr);
    deadcroots = (s->croot * mu_lcr);
    crootX = crootgrowth - deadcroots;

    branchgrowth = (NPP * mu_ab);
    deadbranches = (s->branch * mu_lb);
    branchX = branchgrowth - deadbranches;

    stemgrowth = (NPP * mu_aw);
    deadstems = (s->stem * mu_lw);
    stemX = stemgrowth - deadstems;

    woodX = branchX + stemX + crootX;

    deadsapwood = (mu_lw + p->sapturnover) * s->sapwood;
    sapwoodX += stemgrowth - deadsapwood;

    leaf_material = deadleaves * (1.0 - mu_fmleaf);
    wood_material = deadbranches + deadstems + deadsapwood;

    surf_metab_litter = deadleaves * mu_fmleaf;
    surf_metab_to_active = s->metabsurf * mu_decayrate1 * 0.45;

    soil_metab_litter = deadroots * mu_fmroot;
    soil_metab_to_active = s->metabsoil * mu_decayrate3 * 0.45;

    surf_struct_litter = leaf_material + wood_material;
    structout_surf = s->structsurf * mu_decayrate0;
    surf_struct_to_slow = structout_surf * p->ligshoot * 0.7;
    surf_struct_to_active = structout_surf * (1.0 - p->ligshoot) * 0.55;

    structout_soil = s->structsoil * mu_decayrate2;
    soil_struct_to_slow = structout_soil * p->ligroot * 0.7;
    soil_struct_to_active = structout_soil * (1.0 - p->ligroot) * 0.45;
    soil_struct_litter = deadroots * (1.0 - mu_fmroot) + deadcroots;

    frac_microb_resp = 0.85 - (0.68 * p->finesoil);

    activeout = s->activesoil * mu_decayrate4;
    c_into_active = surf_struct_to_active + soil_struct_to_active + \
                    surf_metab_to_active + soil_metab_to_active + \
                    slow_to_active + passive_to_active;
    active_to_slow = activeout * (1.0 - frac_microb_resp - 0.004);
    active_to_passive = activeout * 0.004;

    slowout = s->slowsoil * mu_decayrate5;
    slow_to_active = slowout * 0.42;
    slow_to_passive = slowout * 0.03;
    c_into_slow = surf_struct_to_slow + soil_struct_to_slow + \
                  active_to_slow;

    passive_to_active = s->passivesoil * mu_decayrate6 * 0.45;
    c_into_passive = active_to_passive + slow_to_passive;

    co2_to_air0 = (structout_surf * \
                   (p->ligshoot * 0.3 + (1.0 - p->ligshoot) * 0.45));
    co2_to_air1 = (structout_soil * \
                  (p->ligroot * 0.3 + (1.0 - p->ligroot) * 0.55));
    co2_to_air2 = s->metabsurf * mu_decayrate1 * 0.55;
    co2_to_air3 = s->metabsoil * mu_decayrate3 * 0.55;
    co2_to_air4 = activeout * frac_microb_resp;
    co2_to_air5 = slowout * 0.55;
    co2_to_air6 = s->passivesoil * mu_decayrate6 * 0.55;

    structsurfX = (surf_struct_litter - \
                  (surf_struct_to_slow + surf_struct_to_active +
                   co2_to_air0));

    structsoilX = (soil_struct_litter - \
                  (soil_struct_to_slow + soil_struct_to_active +
                   co2_to_air1));

    metabsurfX = (surf_metab_litter - \
                  (surf_metab_to_active + co2_to_air2));

    metabsoilX = (soil_metab_litter - \
                  (soil_metab_to_active + co2_to_air3));

    activesoilX = c_into_active - \
                  (active_to_slow + active_to_passive + co2_to_air4);

    slowsoilX = c_into_slow - \
                (slow_to_active + slow_to_passive + co2_to_air5);

    passivesoilX = c_into_passive - \
                    (passive_to_active + co2_to_air6);

    // Update the state
    s->shoot += shootX;
    s->root += rootX;
    s->croot += crootX;
    s->branch += branchX;
    s->stem += stemX;
    s->sapwood += sapwoodX;
    s->metabsoil += metabsoilX;
    s->metabsurf += metabsurfX;
    s->structsoil += structsoilX;
    s->structsurf += structsurfX;
    s->activesoil += activesoilX;
    s->slowsoil += slowsoilX;
    s->passivesoil += passivesoilX;

    // Now solve the N pools using the average NC ratio
    s->shootn = s->shoot * fs->shoot_nc;
    s->rootn = s->root * fs->root_nc;
    s->crootn = s->croot * fs->croot_nc;
    s->branchn = s->branch * fs->branch_nc;
    s->stemn = s->stem * fs->stem_nc;
    s->stemnimm = s->stem * fs->stemnimm_ratio;
    s->stemnmob = s->stem * fs->stemnmob_ratio;
    s->metabsoiln = s->metabsoil * fs->metablsoil_nc;
    s->metabsurfn = s->metabsurf * fs->metabsurf_nc;
    s->structsoiln = s->structsoil * fs->structsoil_nc;
    s->structsurfn = s->structsurf * fs->structsurf_nc;
    s->activesoiln = s->activesoil * fs->activesoil_nc;
    s->slowsoiln = s->slowsoil * fs->slowsoil_nc;
    s->passivesoiln = s->passivesoil * fs->passivesoil_nc;

    // Step 4:  Keep spinning until the slowest C pool (passive) hit equilibrium
    while (TRUE) {
        if (fabs(s->passivesoil - prev_passivec) < 0.05) {
            break;
        } else {
            prev_passivec = s->passivesoil;
            run_sim(cw, c, f, fs, ma, m, p, s, nr);
        }
    }

    fprintf(stderr,
      "Spunup: Plant C - %f, Soil C - %f\n", s->plantc, s->soilc);

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


void zero_fast_spinup_stuff(fast_spinup *fs) {

    int i;

    fs->ndays = 0;
    fs->npp_ss = 0.0;

    for (i = 0; i < 5; i++) {
        fs->alloc[i] = 0.0;
        fs->loss[i] = 0.0;
    }

    fs->alloc[S1] = 0.0;
    fs->alloc[S2] = 0.0;

    for (i = 0; i <= 6; i++) {
        fs->dr[i] = 0.0;
    }

    fs->shoot_nc = 0.0;
    fs->root_nc = 0.0;
    fs->branch_nc = 0.0;
    fs->croot_nc = 0.0;
    fs->stem_nc = 0.0;
    fs->stemnmob_ratio = 0.0;
    fs->stemnimm_ratio = 0.0;
    fs->metablsoil_nc = 0.0;
    fs->metabsurf_nc = 0.0;
    fs->structsoil_nc = 0.0;
    fs->structsurf_nc = 0.0;
    fs->activesoil_nc = 0.0;
    fs->slowsoil_nc = 0.0;
    fs->passivesoil_nc = 0.0;

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

void allocate_numerical_libs_stuff(nrutil *nr) {

    nr->xp = dvector(1, nr->kmax);
    nr->yp = dmatrix(1, nr->N, 1, nr->kmax);
    nr->yscal = dvector(1, nr->N);
    nr->y = dvector(1, nr->N);
    nr->dydx = dvector(1, nr->N);
    nr->ystart = dvector(1, nr->N);
    nr->ak2 = dvector(1, nr->N);
    nr->ak3 = dvector(1, nr->N);
    nr->ak4 = dvector(1, nr->N);
    nr->ak5 = dvector(1, nr->N);
    nr->ak6 = dvector(1, nr->N);
    nr->ytemp = dvector(1, nr->N);
    nr->yerr = dvector(1, nr->N);

    return;
}


void fill_up_solar_arrays(canopy_wk *cw, control *c, met_arrays *ma, params *p) {

    // This is a suprisingly big time hog. So I'm going to unpack it once into
    // an array which we can then access during spinup to save processing time

    int    nyr, doy, hod;
    long   ntimesteps = c->total_num_days * 48;
    double year, sw_rad;

    cw->cz_store = malloc(ntimesteps * sizeof(double));
    if (cw->cz_store == NULL) {
        fprintf(stderr, "malloc failed allocating cz store\n");
        exit(EXIT_FAILURE);
    }

    cw->ele_store = malloc(ntimesteps * sizeof(double));
    if (cw->ele_store == NULL) {
        fprintf(stderr, "malloc failed allocating ele store\n");
        exit(EXIT_FAILURE);
    }

    cw->df_store = malloc(ntimesteps * sizeof(double));
    if (cw->df_store == NULL) {
        fprintf(stderr, "malloc failed allocating df store\n");
        exit(EXIT_FAILURE);
    }

    c->hour_idx = 0;
    for (nyr = 0; nyr < c->num_years; nyr++) {
        year = ma->year[c->hour_idx];
        if (is_leap_year(year))
            c->num_days = 366;
        else
            c->num_days = 365;
        for (doy = 0; doy < c->num_days; doy++) {
            for (hod = 0; hod < c->num_hlf_hrs; hod++) {
                calculate_solar_geometry(cw, p, doy, hod);
                sw_rad = ma->par[c->hour_idx] * PAR_2_SW; /* W m-2 */
                get_diffuse_frac(cw, doy, sw_rad);
                cw->cz_store[c->hour_idx] = cw->cos_zenith;
                cw->ele_store[c->hour_idx] = cw->elevation;
                cw->df_store[c->hour_idx] = cw->diffuse_frac;
                c->hour_idx++;
            }
        }
    }
    return;

}
