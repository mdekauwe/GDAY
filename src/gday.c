/* ============================================================================
* Generic Decomposition And Yield (GDAY) model.
*
* G'DAY is a process based model, which runs on a daily timestep and
* simulates carbon, nutrient and water state and fluxes. See below for model
* description.
*
* Paramaeter descriptions are in gday.h
*
* NOTES:
*   I'm essentially transfering the python to C here...
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   27.02.2015
*
* =========================================================================== */

#include "gday.h"
#include "constants.h"
#include "utilities.h"
#include "litter_production.h"
#include "plant_growth.h"
#include "mate.h"
#include "water_balance.h"
#include "simple_moving_average.h"
#include "soils.h"
#include "optimal_root_model.h"
#include "initialise_model.h"
#include "write_output_file.h"
#include "version.h"
#include "phenology.h"
#include "read_param_file.h"
#include "read_met_file.h"

int main(int argc, char **argv)
{
    int error = 0;

    /*
    ** Setup structures, initialise stuff, e.g. zero fluxes.
    */
    control *c;
    fluxes *f;
    met *m;
    params *p;
    state *s;

    c = (control *)malloc(sizeof (control));
    if (c == NULL) {
        fprintf(stderr, "control structure: Not allocated enough memory!\n");
    	exit(1);
    }

    f = (fluxes *)malloc(sizeof (fluxes));
    if (f == NULL) {
    	fprintf(stderr, "fluxes structure: Not allocated enough memory!\n");
    	exit(1);
    }

    m = (met *)malloc(sizeof (met));
    if (m == NULL) {
    	fprintf(stderr, "met structure: Not allocated enough memory!\n");
    	exit(1);
    }

    p = (params *)malloc(sizeof (params));
    if (p == NULL) {
    	fprintf(stderr, "params structure: Not allocated enough memory!\n");
    	exit(1);
    }

    s = (state *)malloc(sizeof (state));
    if (s == NULL) {
    	fprintf(stderr, "state structure: Not allocated enough memory!\n");
    	exit(1);
    }

    initialise_control(c);
    initialise_params(p);
    initialise_fluxes(f);
    initialise_state(s);

    clparser(argc, argv, c);
    /*
    ** Read .ini parameter file and meterological data
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

    read_met_data(argv, c, m);

    if (c->spin_up)
        spin_up_pools(c, f, m, p, s);
    else
        run_sim(c, f, m, p, s);

    /* clean up */
    fclose(c->ofp);
    fclose(c->ifp);
    if (c->output_ascii == FALSE) {
        fclose(c->ofp_hdr);
    }

    free(c);
    free(f);
    free(m->year);
    free(m->prjday);
    free(m->sw_rad);
    free(m->tair);
    free(m->rain);
    free(m->tsoil);
    free(m->tam);
    free(m->tpm);
    free(m->vpd_am);
    free(m->vpd_pm);
    free(m->vpd_avg);
    free(m->co2);
    free(m->ndep);
    free(m->wind);
    free(m->press);
    free(m->par);
    free(m->wind_am);
    free(m->wind_pm);
    free(m->sw_rad_am);
    free(m->sw_rad_pm);
    free(m);
    free(p);
    free(s);

    exit(EXIT_SUCCESS);
}



void run_sim(control *c, fluxes *f, met *m, params *p, state *s){

    int nyr, doy, window_size, i;
    int project_day = 0;

    double fdecay, rdecay, current_limitation, nitfac, year;

    /* potentially allocating 1 extra spot, but will be fine as we always
       index by num_days */
    double *day_length = NULL;
    if ((day_length = (double *)calloc(366, sizeof(double))) == NULL) {
        fprintf(stderr,"Error allocating space for day_length\n");
		exit(EXIT_FAILURE);
    }

    if (c->deciduous_model) {
        if (s->max_lai < -900.){
            /* initialise to something really low */
            s->max_lai = 0.01;
            s->max_shoot = 0.01;
        }

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
    if (c->print_options == DAILY && c->spin_up == FALSE) {
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


    /* Set up SMA
    **  - If we don't have any information about the N & water limitation, i.e.
    **    as would be the case with spin-up, assume that there is no limitation
    **    to begin with.
    */
    if (s->prev_sma < -900)
        s->prev_sma = 1.0;
    if (s->grw_seas_stress < -900)
        s->grw_seas_stress = 1.0;

    /*
        Window size = root lifespan in days...
        For deciduous species window size is set as the length of the
        growing season in the main part of the code
    */
    window_size = (int)(1.0 / p->rdecay * NDAYS_IN_YR);
    sma_obj *hw = sma(SMA_NEW, window_size).handle;
    for (i = 0; i < window_size; i++) {
        sma(SMA_ADD, hw, s->prev_sma);
    }

    /*
        params are defined in per year, needs to be per day. Important this is
        done here as rate constants elsewhere in the code are assumed to be in
        units of days not years
    */
    correct_rate_constants(p, FALSE);
    day_end_calculations(c, p, s, -99, TRUE);

    initialise_soil_moisture_parameters(c, p);
    s->pawater_root = p->wcapac_root;
    s->pawater_topsoil = p->wcapac_topsoil;
    
    s->lai = MAX(0.01, (p->sla * M2_AS_HA / KG_AS_TONNES /
                        p->cfracts * s->shoot));

    /* ====================== **
    **   Y E A R    L O O P   **
    ** ====================== */
    project_day = 0;
    for (nyr = 0; nyr < c->num_years; nyr++) {
        year = m->year[project_day];
        if (is_leap_year(year))
            c->num_days = 366;
        else
            c->num_days = 365;

        calculate_daylength(c->num_days, p->latitude, *(&day_length));

        if (c->deciduous_model) {
            phenology(c, f, m, p, s, day_length, project_day);

            /* Change window size to length of growing season */
            sma(SMA_FREE, hw);
            hw = sma(SMA_NEW, p->growing_seas_len).handle;
            for (i = 0; i < p->growing_seas_len; i++) {
                sma(SMA_ADD, hw, 1.0); /* don't rely on previous year */
            }

            zero_stuff(c, s);
        }
        /* =================== **
        **   D A Y   L O O P   **
        ** =================== */

        for (doy = 0; doy < c->num_days; doy++) {

            calculate_litterfall(c, f, p, s, doy, &fdecay, &rdecay);

            calc_day_growth(c, f, m, p, s, project_day, day_length[doy],
                            doy, fdecay, rdecay);

            calculate_csoil_flows(c, f, p, s, m->tsoil[project_day]);
            calculate_nsoil_flows(c, f, p, s, m->ndep[project_day]);

            /*printf("%f\n", f->gpp*100.);*/

            /* update stress SMA */
            if (c->deciduous_model && s->leaf_out_days[doy] > 0.0) {
                 /*Allocation is annually for deciduous "tree" model, but we
                   need to keep a check on stresses during the growing season
                   and the LAI figure out limitations during leaf growth period.
                   This also applies for deciduous grasses, need to do the
                   growth stress calc for grasses here too. */
                current_limitation = calculate_growth_stress_limitation(p, s);
                sma(SMA_ADD, hw, current_limitation);
                s->prev_sma = sma(SMA_MEAN, hw).sma;
            } else if (c->deciduous_model == FALSE) {
                current_limitation = calculate_growth_stress_limitation(p, s);
                sma(SMA_ADD, hw, current_limitation);
                s->prev_sma = sma(SMA_MEAN, hw).sma;
            }

            /* if grazing took place need to reset "stress" running mean calculation
               for grasses */
            if (c->grazing == 2) {
                sma(SMA_FREE, hw);
                hw = sma(SMA_NEW, window_size).handle;
                for (i = 0; i < window_size; i++) {
                    sma(SMA_ADD, hw, s->prev_sma);
                }
            }


            /* Turn off all N calculations */
            if (c->ncycle == FALSE)
                reset_all_n_pools_and_fluxes(f, s);

            /* calculate C:N ratios and increment annual flux sum */
            day_end_calculations(c, p, s, c->num_days, FALSE);

            if (c->print_options == DAILY && c->spin_up == FALSE) {
                if(c->output_ascii)
                    write_daily_outputs_ascii(c, f, s, year, doy+1);
                else
                    write_daily_outputs_binary(c, f, s, year, doy+1);
            }


            /* check the daily water balance */
            /*check_water_balance(project_day); */

            project_day++;
            /* ======================= **
            **   E N D   O F   D A Y   **
            ** ======================= */
        }

        /* Allocate stored C&N for the following year */
        if (c->deciduous_model) {
            calculate_average_alloc_fractions(f, s, p->growing_seas_len);
            allocate_stored_c_and_n(f, p, s);
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

    return;


}

void spin_up_pools(control *c, fluxes *f, met *m, params *p, state *s){
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
    int i;
    /* check for convergences in units of kg/m2 */
    double conv = TONNES_HA_2_KG_M2;


    /* Final state + param file */
    open_output_file(c, c->out_param_fname, &(c->ofp));

    fprintf(stderr, "Spinning up the model...\n");
    while (TRUE) {
        if (fabs((prev_plantc*conv) - (s->plantc*conv)) < tol &&
            fabs((prev_soilc*conv) - (s->soilc*conv)) < tol) {
            break;
        } else {
            prev_plantc = s->plantc;
            prev_soilc = s->soilc;

            /* 1000 years (50 yrs x 20 cycles) */
            for (i = 0; i < 20; i++)
                run_sim(c, f, m, p, s); /* run GDAY */

            /* Have we reached a steady state? */
            fprintf(stderr, "Spinup: Plant C - %f, Soil C - %f\n",
                    s->plantc, s->soilc);
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
    fprintf(stderr, "\t%s [options]\n", argv[0]);
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
    /* If the N-Cycle is turned off the way I am implementing this is to
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
    s->max_lai = 0.0;
    s->max_shoot = 0.0;
    s->cstore = 0.0;
    s->nstore = 0.0;
    s->anpp = 0.0;
    s->grw_seas_stress = 1.0;

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
