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

    // zeroed: the canopy air iteration starts from state kept in cw
    cw = (canopy_wk *)calloc(1, sizeof(canopy_wk));
    if (cw == NULL) {
        fprintf(stderr, "canopy wk structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    f = (fluxes *)malloc(sizeof(fluxes));
    if (f == NULL) {
    	fprintf(stderr, "fluxes structure: Not allocated enough memory!\n");
    	exit(EXIT_FAILURE);
    }

    ma = (met_arrays *)calloc(1, sizeof(met_arrays));
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
    if (fs == NULL) {
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
    zero_fast_spinup_stuff(fs);

    clparser(argc, argv, c);
    strcpy(c->git_code_ver, build_git_sha);
    if (c->PRINT_GIT) {
        fprintf(stderr, "\n%s\n", c->git_code_ver);
        exit(EXIT_SUCCESS);
    }

    /*
     * Read .ini parameter file and meterological data
     */
    error = parse_ini_file(c, p, s);
    if (error != 0) {
        prog_error("Error reading .INI file on line", __LINE__);
    }

    /* House keeping! */
    if (c->water_balance == HYDRAULICS && c->sub_daily == FALSE) {
        fprintf(stderr, "You can't run the hydraulics model with daily flag\n");
        exit(EXIT_FAILURE);
    }
    if (c->canopy_evap_model == CANOPY_EVAP_JULES && !c->sub_daily) {
        fprintf(stderr, "canopy_evap_model = jules needs sub_daily = true\n");
        exit(EXIT_FAILURE);
    }
    if (c->soil_evap_model == SOIL_EVAP_OR && c->water_balance != HYDRAULICS) {
        /* needs the conductivity of the top SPA layer */
        fprintf(stderr, "soil_evap_model = or needs water_balance = hydraulics\n");
        exit(EXIT_FAILURE);
    }
    if (c->pcycle && c->deciduous_model) {
        fprintf(stderr, "pcycle isn't implemented for the deciduous model "
                "yet\n");
        exit(EXIT_FAILURE);
    }
    if (c->p_limit_photo && (c->modeljm < 1 || c->modeljm > 2 ||
                             c->ps_pathway != C3)) {
        fprintf(stderr, "p_limit_photo caps the N-based Vcmax/Jmax, so needs "
                "modeljm = 1 or 2 and C3 photosynthesis\n");
        exit(EXIT_FAILURE);
    }
    if (c->water_store) {
        fprintf(stderr, "water_store (plant capacitance) was part of the Emax "
                "scheme, which gs_opt replaced (git tag last-emax)\n");
        exit(EXIT_FAILURE);
    }
    if (c->water_balance == HYDRAULICS && !c->sub_daily) {
        fprintf(stderr, "water_balance = hydraulics (gs_opt) needs "
                "sub_daily = true\n");
        exit(EXIT_FAILURE);
    }
    if (c->water_balance == HYDRAULICS &&
        (p->p50 >= 0.0 || p->p88 >= p->p50)) {
        fprintf(stderr, "gs_opt needs p88 < p50 < 0 (MPa)\n");
        exit(EXIT_FAILURE);
    }
    if (c->plant_segments == N_PLANT_SEG) {
        /* each segment's effective traits (unset = the whole plant's) */
        double fr[N_PLANT_SEG] = {p->seg_frac_root, p->seg_frac_stem,
                                  p->seg_frac_leaf};
        double s50[N_PLANT_SEG] = {p->p50_root, p->p50_stem, p->p50_leaf};
        double s88[N_PLANT_SEG] = {p->p88_root, p->p88_stem, p->p88_leaf};
        const char *nm[N_PLANT_SEG] = {"root", "stem", "leaf"};
        int k;
        for (k = 0; k < N_PLANT_SEG; k++) {
            double a50 = s50[k] < -900.0 ? p->p50 : s50[k];
            double a88 = s88[k] < -900.0 ? p->p88 : s88[k];
            if (fr[k] <= 0.0 || a50 >= 0.0 || a88 >= a50) {
                fprintf(stderr, "plant_segments: the %s segment needs "
                        "seg_frac > 0 and p88 < p50 < 0 (MPa), got "
                        "seg_frac %g, p50 %g, p88 %g\n", nm[k], fr[k], a50,
                        a88);
                exit(EXIT_FAILURE);
            }
        }
    }

    if (c->gs_model == GS_OPT) {
        /* daily (MATE) profit maximisation on the bucket */
        if (c->sub_daily || c->water_balance != BUCKET || c->ps_pathway != C3) {
            fprintf(stderr, "gs_model = gs_opt is the daily C3 (MATE) model "
                    "with water_balance = bucket; the sub-daily model uses "
                    "gs_opt through water_balance = hydraulics\n");
            exit(EXIT_FAILURE);
        }
        if (c->nonstomatal_limitation) {
            fprintf(stderr, "gs_opt has no non-stomatal limitation, set "
                    "nonstomatal_limitation = false\n");
            exit(EXIT_FAILURE);
        }
        if (p->p50 >= 0.0 || p->p88 >= p->p50 || p->kp <= 0.0) {
            fprintf(stderr, "gs_opt needs kp > 0 and p88 < p50 < 0 (MPa)\n");
            exit(EXIT_FAILURE);
        }
        /* bucket water potential from the retention curve */
        setup_soil_retention(c, p);
    }

    if (c->water_balance == HYDRAULICS) {
        allocate_numerical_libs_stuff(nr);
        initialise_roots(f, p, s);
        setup_hydraulics_arrays(f, p, s);
    }

    /* PLUMBER2 style netCDF (as JULES reads) or GDAY's ascii format */
    if (is_netcdf_file(c->met_fname)) {
        read_met_data_netcdf(argv, c, ma, p);
    } else if (c->sub_daily) {
        read_subdaily_met_data(argv, c, ma);
    } else {
        read_daily_met_data(argv, c, ma);
    }
    if (c->sub_daily) {
        fill_up_solar_arrays(cw, c, ma, p);
    }


    if (c->spin_up) {
        spin_up_pools(cw, c, f, fs, ma, m, p, s, nr);
    } else {
        run_sim(cw, c, f, fs, ma, m, p, s, nr);
    }

    if (c->water_balance == HYDRAULICS && c->n_supply_bind > 0) {
        fprintf(stderr, "gs_opt: the soil supply limited %ld leaf steps\n",
                c->n_supply_bind);
    }

    /* clean up - not every file is opened in every run mode */
    if (c->ofp != NULL) {
        fclose(c->ofp);
    }
    if (c->ofp_sd != NULL) {
        fclose(c->ofp_sd);
    }
    if (c->ifp != NULL) {
        fclose(c->ifp);
    }
    if (c->ofp_hdr != NULL) {
        fclose(c->ofp_hdr);
    }

    free(ma->year);
    free(ma->tair);
    free(ma->rain);
    free(ma->tsoil);
    free(ma->lwdown);
    free(ma->co2);
    free(ma->ndep);
    free(ma->nfix);
    free(ma->lai);
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
        free(ma->tday);
        free(ma->vpd_am);
        free(ma->vpd_pm);
        free(ma->wind_am);
        free(ma->wind_pm);
        free(ma->par_am);
        free(ma->par_pm);
    }
    free(s->day_length);
    free(cw);
    free(c);
    free(ma);
    free(m);
    free(p);
    free(s);
    free(f);
    free(fs);
    free(nr);

    exit(EXIT_SUCCESS);
}

double total_ecosystem_n(state *s) {
    /* all the N in the plant, its store, litter, SOM and inorganic pools
       (t ha-1) */
    return (s->shootn + s->rootn + s->crootn + s->branchn + s->stemnimm +
            s->stemnmob + s->nstore + s->structsurfn + s->structsoiln +
            s->metabsurfn + s->metabsoiln + s->activesoiln + s->slowsoiln +
            s->passivesoiln + s->inorgn);
}

void check_n_balance(fluxes *f, state *s, double n_start, int year,
                     int doy) {
    /* the day's change in ecosystem N must equal inputs - losses */
    double err = total_ecosystem_n(s) - n_start - (f->ninflow - f->nloss);

    if (fabs(err) > 1E-9 * MAX(1.0, n_start)) {
        fprintf(stderr, "N balance not closed: %d doy %d error %.3e t ha-1\n",
                year, doy + 1, err);
        exit(EXIT_FAILURE);
    }
}

/*
** During spin-up the mineral P pools are held at the site's values (with
** deposition and little loss they have no steady state and would drift up
** over the millennia of a spin-up), so the plant and organic P come to
** equilibrium with the site's mineral P
*/
static int    hold_mineral_p = FALSE;
static double held_mineral_p[5];

static void hold_mineral_p_pools(params *p, state *s) {
    /* labile & sorbed put on the Langmuir isotherm from the labile P */
    s->inorgsorbp = p->smax * s->inorglabp / (p->ks + s->inorglabp);
    held_mineral_p[0] = s->inorglabp;
    held_mineral_p[1] = s->inorgsorbp;
    held_mineral_p[2] = s->inorgssorbp;
    held_mineral_p[3] = s->inorgoccp;
    held_mineral_p[4] = s->inorgparp;
    hold_mineral_p = TRUE;
}

static void restore_mineral_p_pools(state *s) {
    s->inorglabp = held_mineral_p[0];
    s->inorgsorbp = held_mineral_p[1];
    s->inorgssorbp = held_mineral_p[2];
    s->inorgoccp = held_mineral_p[3];
    s->inorgparp = held_mineral_p[4];
}

double total_ecosystem_p(state *s) {
    /* all the P in the plant, litter, SOM and mineral pools, including the
       parent material (t ha-1) */
    return (s->shootp + s->rootp + s->crootp + s->branchp + s->stempimm +
            s->stempmob + s->pstore + s->structsurfp + s->structsoilp +
            s->metabsurfp + s->metabsoilp + s->activesoilp + s->slowsoilp +
            s->passivesoilp + s->inorglabp + s->inorgsorbp + s->inorgssorbp +
            s->inorgoccp + s->inorgparp);
}

void check_p_balance(fluxes *f, state *s, double p_start, int year,
                     int doy) {
    /* the day's change in ecosystem P must equal deposition - leaching */
    double err = total_ecosystem_p(s) - p_start - (f->p_atm_dep - f->ploss);

    if (fabs(err) > 1E-9 * MAX(1.0, p_start)) {
        fprintf(stderr, "P balance not closed: %d doy %d error %.3e t ha-1\n",
                year, doy + 1, err);
        exit(EXIT_FAILURE);
    }
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

    if (c->fixed_lai) {
        s->lai = p->fix_lai;
    } else if (c->prescribed_lai) {
        /* reset each day in the day loop */
        s->lai = ma->lai[0];
    } else {
        s->lai = MAX(0.01, (p->sla * M2_AS_HA / KG_AS_TONNES /
                            p->cfracts * s->shoot));
    }

    if (c->water_balance == HYDRAULICS) {
        double root_zone_total, water_content;

        // Update the soil water storage
        root_zone_total = 0.0;
        for (i = 0; i < p->soil_layers; i++) {

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
            }

            /* Hurricane? NB. hurricane_doy is 1-based like disturbance_doy */
            if (c->hurricane &&
                p->hurricane_yr == year &&
                p->hurricane_doy == doy+1) {
                hurricane(f, p, s);
            }


            if (c->prescribed_lai) {
                set_prescribed_lai(c, ma, s);
            }

#ifdef CHECK_NUTRIENT_BALANCE
            double n_start = total_ecosystem_n(s);
            double p_start = total_ecosystem_p(s);
#endif
            calc_day_growth(cw, c, f, fs, ma, m, nr, p, s, s->day_length[doy],
                            doy, fdecay, rdecay);

            //printf("%d %f %f\n", doy, f->gpp*100, s->lai);
            calculate_csoil_flows(c, f, fs, p, s, m->tsoil, doy);
            calculate_nsoil_flows(c, f, p, s, doy);
            if (c->pcycle) {
                calculate_psoil_flows(c, f, p, s, doy);
            }
#ifdef CHECK_NUTRIENT_BALANCE
            if (c->ncycle) {
                check_n_balance(f, s, n_start, year, doy);
            }
            if (c->pcycle) {
                check_p_balance(f, s, p_start, year, doy);
            }
#endif
            if (hold_mineral_p) {
                restore_mineral_p_pools(s);
            }

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
                write_daily_outputs_ascii(c, cw, f, p, s, year, doy+1);
            } else if (c->print_options == DAILY && c->spin_up == FALSE) {
                if(c->output_ascii)
                    write_daily_outputs_ascii(c, cw, f, p, s, year, doy+1);
                else
                    write_daily_outputs_binary(c, f, s, year, doy+1);
            }

            // Store the time-varying variables for the SAS spin-up
            if (c->spinup_method == SAS) {
                accumulate_fast_spinup_stuff(fs, f, s);
            }
            c->day_idx++;
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
    if (c->pcycle) {
        hold_mineral_p_pools(p, s);
    }

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

    hold_mineral_p = FALSE;
    write_final_state(c, p, s);

    return;
}

void sas_spinup(canopy_wk *cw, control *c, fluxes *f, fast_spinup *fs,
                met_arrays *ma, met *m, params *p, state *s, nrutil *nr) {
    //
    // Semi-analytical solution (SAS) to accelerate model spin-up of
    // carbon–nitrogen pools, following Xia et al. (2013) GMD.
    //
    // 1. Run the model until the fast plant pools (leaves, fine roots) and
    //    hence NPP and allocation are ~stable.
    // 2. Set the woody, litter and SOM pools to the steady state of their
    //    (linear) dynamics, dX/dt = inputs - k X = 0, using the mean inputs
    //    and turnover rates over the last cycle.
    // 3. Run the model until plant and soil C are stable, which also lets
    //    the N pools and any feedbacks (e.g. allometry) settle.
    //
    // NB. root exudation inputs to the active pool are not included in the
    //     analytical step, step 3 takes care of that. Pools with a zero
    //     turnover rate (e.g. crdecay = 0) have no steady state and are left
    //     as simulated.
    //

    double cleaf0, croot0, rel_change, ratio;
    double ndays, k[7], in_ss, in_ms, in_st, in_mt;
    double bd, wd, crd, sapd, cpbranch, cpstem, cpcroot;
    double to_active, to_slow, frac_microb_resp, x, a, b, cp;
    double prev_plantc, prev_soilc;
    int    i, iter;

    // Steps 1 & 2 are repeated, as NPP responds (via N mineralisation) to
    // the soil pools just solved for, until the analytical solution settles
    run_sim(cw, c, f, fs, ma, m, p, s, nr);
    prev_plantc = -999.9;
    prev_soilc = -999.9;
    for (iter = 0; iter < 50; iter++) {

    // Step 1: fast plant pools to steady state (summed relative change of
    //         leaf and fine root C over a cycle < 1%)
    while (TRUE) {
        cleaf0 = s->shoot;
        croot0 = s->root;

        zero_fast_spinup_stuff(fs);
        run_sim(cw, c, f, fs, ma, m, p, s, nr);

        rel_change = (fabs(s->shoot - cleaf0) / MAX(s->shoot, 1E-06) +
                      fabs(s->root - croot0) / MAX(s->root, 1E-06));
        if (rel_change < 0.01) {
            break;
        }
    }

    // Step 2: analytical steady state from the means of the last cycle
    ndays = (double)fs->ndays;
    for (i = 0; i < 7; i++) {
        k[i] = fs->dr[i] / ndays;
    }
    cpbranch = fs->cpbranch / ndays;
    cpstem = fs->cpstem / ndays;
    cpcroot = fs->cpcroot / ndays;

    // run_sim leaves the rate constants in per year units
    bd = p->bdecay / NDAYS_IN_YR;
    wd = p->wdecay / NDAYS_IN_YR;
    crd = p->crdecay / NDAYS_IN_YR;
    sapd = p->sapturnover / NDAYS_IN_YR;

    // litter inputs; at steady state woody litter = woody growth
    in_ss = fs->surf_struct_litter / ndays;
    in_ms = fs->surf_metab_litter / ndays;
    in_st = fs->soil_struct_litter / ndays;
    in_mt = fs->soil_metab_litter / ndays;

    // woody pools, N scaled so their N:C is unchanged
    if (bd > 0.0 && s->branch > 0.0) {
        ratio = (cpbranch / bd) / s->branch;
        s->branch *= ratio;
        s->branchn *= ratio;
        in_ss += cpbranch - fs->deadbranch / ndays;
    }
    if (wd > 0.0 && s->stem > 0.0) {
        ratio = (cpstem / wd) / s->stem;
        s->stem *= ratio;
        s->stemnimm *= ratio;
        s->stemnmob *= ratio;
        s->stemn = s->stemnimm + s->stemnmob;
        s->sapwood = cpstem / (wd + sapd);
        in_ss += cpstem - fs->deadstems / ndays;
    }
    if (crd > 0.0 && s->croot > 0.0) {
        ratio = (cpcroot / crd) / s->croot;
        s->croot *= ratio;
        s->crootn *= ratio;
        in_st += cpcroot - fs->deadcroots / ndays;
    }

    // litter pools, all their outflow is k X = inputs
    s->structsurf = in_ss / k[0];
    s->metabsurf = in_ms / k[1];
    s->structsoil = in_st / k[2];
    s->metabsoil = in_mt / k[3];

    // litter C passed on to the SOM pools (same fractions as soils.c)
    to_active = (in_ss * (1.0 - p->ligshoot) * 0.55 +
                 in_st * (1.0 - p->ligroot) * 0.45 +
                 in_ms * 0.45 + in_mt * 0.45);
    to_slow = in_ss * p->ligshoot * 0.7 + in_st * p->ligroot * 0.7;

    // SOM pools are coupled. Writing a, b, cp for the outflows kX of the
    // active, slow and passive pools:
    //   a  = to_active + 0.42 b + 0.45 cp
    //   b  = to_slow + x a,          x = 1 - frac_microb_resp - 0.004
    //   cp = 0.004 a + 0.03 b
    frac_microb_resp = 0.85 - (0.68 * p->finesoil);
    x = 1.0 - frac_microb_resp - 0.004;
    if (c->passiveconst) {
        // passive pool is fixed, so its outflow is known
        cp = k[6] * s->passivesoil;
        a = (to_active + 0.42 * to_slow + 0.45 * cp) / (1.0 - 0.42 * x);
        b = to_slow + x * a;
    } else {
        a = ((to_active + (0.42 + 0.45 * 0.03) * to_slow) /
             (1.0 - 0.42 * x - 0.45 * (0.004 + 0.03 * x)));
        b = to_slow + x * a;
        cp = 0.004 * a + 0.03 * b;
        s->passivesoil = cp / k[6];
    }
    s->activesoil = a / k[4];
    s->slowsoil = b / k[5];

    // N pools from the mean N:C ratios of the last cycle
    s->metabsoiln = s->metabsoil * fs->metablsoil_nc / ndays;
    s->metabsurfn = s->metabsurf * fs->metabsurf_nc / ndays;
    s->structsoiln = s->structsoil * fs->structsoil_nc / ndays;
    s->structsurfn = s->structsurf * fs->structsurf_nc / ndays;
    s->activesoiln = s->activesoil * fs->activesoil_nc / ndays;
    s->slowsoiln = s->slowsoil * fs->slowsoil_nc / ndays;
    if (c->passiveconst == FALSE) {
        s->passivesoiln = s->passivesoil * fs->passivesoil_nc / ndays;
    }

    // derived totals
    day_end_calculations(c, p, s, -99, TRUE);
    fprintf(stderr, "SAS spinup (analytical): Plant C - %f, Soil C - %f\n",
            s->plantc, s->soilc);

    if (fabs(s->plantc - prev_plantc) < 1E-03 * s->plantc &&
        fabs(s->soilc - prev_soilc) < 1E-03 * s->soilc) {
        break;
    }
    prev_plantc = s->plantc;
    prev_soilc = s->soilc;

    } /* end of analytical iterations */

    // Step 3: run until plant & soil C change by < 0.01% per cycle
    while (TRUE) {
        prev_plantc = s->plantc;
        prev_soilc = s->soilc;
        zero_fast_spinup_stuff(fs);
        run_sim(cw, c, f, fs, ma, m, p, s, nr);
        fprintf(stderr, "SAS spinup: Plant C - %f, Soil C - %f\n",
                s->plantc, s->soilc);

        if (fabs(s->plantc - prev_plantc) < 1E-04 * s->plantc &&
            fabs(s->soilc - prev_soilc) < 1E-04 * s->soilc) {
            break;
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
    for (i = 0; i < 7; i++) {
        fs->dr[i] = 0.0;
    }
    fs->surf_struct_litter = 0.0;
    fs->surf_metab_litter = 0.0;
    fs->soil_struct_litter = 0.0;
    fs->soil_metab_litter = 0.0;
    fs->cpbranch = 0.0;
    fs->cpstem = 0.0;
    fs->cpcroot = 0.0;
    fs->deadbranch = 0.0;
    fs->deadstems = 0.0;
    fs->deadcroots = 0.0;
    fs->metablsoil_nc = 0.0;
    fs->metabsurf_nc = 0.0;
    fs->structsoil_nc = 0.0;
    fs->structsurf_nc = 0.0;
    fs->activesoil_nc = 0.0;
    fs->slowsoil_nc = 0.0;
    fs->passivesoil_nc = 0.0;

    return;
}

void accumulate_fast_spinup_stuff(fast_spinup *fs, fluxes *f, state *s) {
    /* daily sums of woody growth/turnover and of the litter & SOM N:C
       ratios (0 when a pool is empty) */

    fs->ndays++;
    fs->cpbranch += f->cpbranch;
    fs->cpstem += f->cpstem;
    fs->cpcroot += f->cpcroot;
    fs->deadbranch += f->deadbranch;
    fs->deadstems += f->deadstems;
    fs->deadcroots += f->deadcroots;
    if (s->metabsoil > 0.0)
        fs->metablsoil_nc += s->metabsoiln / s->metabsoil;
    if (s->metabsurf > 0.0)
        fs->metabsurf_nc += s->metabsurfn / s->metabsurf;
    if (s->structsoil > 0.0)
        fs->structsoil_nc += s->structsoiln / s->structsoil;
    if (s->structsurf > 0.0)
        fs->structsurf_nc += s->structsurfn / s->structsurf;
    if (s->activesoil > 0.0)
        fs->activesoil_nc += s->activesoiln / s->activesoil;
    if (s->slowsoil > 0.0)
        fs->slowsoil_nc += s->slowsoiln / s->slowsoil;
    if (s->passivesoil > 0.0)
        fs->passivesoil_nc += s->passivesoiln / s->passivesoil;

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
        p->prateuptake *= NDAYS_IN_YR;
        p->prateloss *= NDAYS_IN_YR;
        p->puptakez *= NDAYS_IN_YR;
        p->p_atm_deposition *= NDAYS_IN_YR;
        p->p_rate_par_weather *= NDAYS_IN_YR;
        p->max_p_biochemical *= NDAYS_IN_YR;
        p->rate_sorb_ssorb *= NDAYS_IN_YR;
        p->rate_ssorb_occ *= NDAYS_IN_YR;
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
        p->prateuptake /= NDAYS_IN_YR;
        p->prateloss /= NDAYS_IN_YR;
        p->puptakez /= NDAYS_IN_YR;
        p->p_atm_deposition /= NDAYS_IN_YR;
        p->p_rate_par_weather /= NDAYS_IN_YR;
        p->max_p_biochemical /= NDAYS_IN_YR;
        p->rate_sorb_ssorb /= NDAYS_IN_YR;
        p->rate_ssorb_occ /= NDAYS_IN_YR;
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

    /* total plant, soil & litter phosphorus */
    if (c->pcycle) {
        s->shootpc = s->shoot > 0.0 ? s->shootp / s->shoot : 0.0;
        s->rootpc = s->root > 0.0 ? MAX(0.0, s->rootp / s->root) : 0.0;
        s->inorgavlp = s->inorglabp + s->inorgsorbp;
        s->inorgp = s->inorgavlp + s->inorgssorbp + s->inorgoccp +
                    s->inorgparp;
        s->soilp = s->inorgp + s->activesoilp + s->slowsoilp +
                   s->passivesoilp;
        s->litterpag = s->structsurfp + s->metabsurfp;
        s->litterpbg = s->structsoilp + s->metabsoilp;
        s->litterp = s->litterpag + s->litterpbg;
        s->stemp = s->stempimm + s->stempmob;
        s->plantp = s->shootp + s->rootp + s->crootp + s->branchp +
                    s->stemp;
        s->totalp = s->plantp + s->pstore + s->litterp + s->soilp;
    }

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
        // wind and VPD floored (0.1 m s-1, 0.05 kPa): flux site forcing can
        // have calm or saturated steps, which would divide by zero in the
        // aerodynamic conductances and the Medlyn gs model
        m->wind = MAX(WIND_MIN, ma->wind[c->hour_idx]);
        m->press = ma->press[c->hour_idx] * KPA_2_PA;
        m->vpd = MAX(VPD_MIN, ma->vpd[c->hour_idx]) * KPA_2_PA;
        m->tair = ma->tair[c->hour_idx];
        m->tsoil = ma->tsoil[c->hour_idx];
        m->lwdown = ma->lwdown != NULL ? ma->lwdown[c->hour_idx] : -999.9;
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
        m->lwdown_am = ma->lwdown_am != NULL ? ma->lwdown_am[c->day_idx]
                                             : -999.9;
        m->lwdown_pm = ma->lwdown_pm != NULL ? ma->lwdown_pm[c->day_idx]
                                             : -999.9;
        m->rain = ma->rain[c->day_idx];
        m->vpd_am = MAX(VPD_MIN, ma->vpd_am[c->day_idx]) * KPA_2_PA;
        m->vpd_pm = MAX(VPD_MIN, ma->vpd_pm[c->day_idx]) * KPA_2_PA;
        m->wind_am = MAX(WIND_MIN, ma->wind_am[c->day_idx]);
        m->wind_pm = MAX(WIND_MIN, ma->wind_pm[c->day_idx]);
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

void set_prescribed_lai(control *c, met_arrays *ma, state *s) {
    /*
        Today's LAI from the met forcing, the daily mean for sub-daily input
        (c->hour_idx points at the first timestep of the day)
    */
    int    i;
    double sum = 0.0;

    if (c->sub_daily) {
        for (i = 0; i < c->num_hlf_hrs; i++) {
            sum += ma->lai[c->hour_idx + i];
        }
        s->lai_prescribed = sum / (double)c->num_hlf_hrs;
    } else {
        s->lai_prescribed = ma->lai[c->day_idx];
    }
    s->lai = s->lai_prescribed;

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
