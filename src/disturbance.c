#include "disturbance.h"


void figure_out_years_with_disturbances(control *c, met_arrays *ma, params *p,
                                        int **yrs, int *cnt) {
    int nyr, year_of_disturbance, yrs_till_event, prjday, year;

    if (p->burn_specific_yr < -900.0) {
        year_of_disturbance = p->burn_specific_yr;
        (*yrs)[0] = p->burn_specific_yr;
    } else {
        yrs_till_event = time_till_next_disturbance();
        year = (int)ma->year[c->day_idx];
        /*year_of_disturbance = year + yrs_till_event; */
        year_of_disturbance = 1996;

        /* figure out the years of the disturbance events  */
        *cnt = 0;
        prjday = 0;

        for (nyr = 0; nyr < c->num_years - 1; nyr++) {
            year = (int)ma->year[c->day_idx];
            if (is_leap_year(year))
                prjday+=366;
            else
                prjday+=365;

            if (year == year_of_disturbance) {
                yrs_till_event = time_till_next_disturbance();
                if (*cnt == 0) {
                    (*yrs)[0] = year_of_disturbance;
                    *cnt += 1;
                } else {
                    *cnt += 1;
                    if ((yrs = (int **)realloc(yrs, (1 + *cnt) * sizeof(int))) == NULL) {
                        fprintf(stderr,"Error resizing years array\n");
                		exit(EXIT_FAILURE);
                    }
                    (*yrs)[*cnt] = year_of_disturbance;
                }

                /* See if there is another event? */
                year_of_disturbance = year + yrs_till_event;
            }
        }
    }

    return;
}

int time_till_next_disturbance(void) {
    /* calculate the number of years until a disturbance event occurs
    assuming a return interval of X years

    - section 3.4.1 D. Knuth, The Art of Computer Programming.

    Parameters
    ----------
    return_interval : int/float
        interval disturbance return at in years
    */
    /* rate = 1.0 / p->return_interval; */

    /*return int(-log(1.0 - random.random()) / rate); */

    return (11);
}

int check_for_fire(control *c, fluxes *f, params *p, state *s, int year,
                   int *distrubance_yrs, int num_disturbance_yrs) {
    /* Check if the current year has a fire, if so "burn" and then
       return an indicator to tell the main code to reset the stress stream */
    int fire_found = FALSE;
    int nyr;

    for (nyr = 0; nyr < num_disturbance_yrs; nyr++) {
        if (year == distrubance_yrs[nyr]) {
            fire_found = TRUE;
        }
    }

    return (fire_found);
}

void fire(control *c, fluxes *f, params *p, state *s) {
    /*
    Fire...

    * 100 percent of aboveground biomass
    * 100 percent of surface litter
    * 50 percent of N volatilized to the atmosphere
    * 50 percent of N returned to inorgn pool
    * 100 percent of N returned to inorgmin pool  Wang 2015 Chemosphere DOI:10.1016/j.chemosphere.2014.01.084
    * Coarse roots are not damaged by fire!

    vaguely following ...
    http://treephys.oxfordjournals.org/content/24/7/765.full.pdf
    */
    double totaln, totalp;

    totaln = s->branchn + s->shootn + s->stemn + s->structsurfn;
    s->inorgn += totaln / 2.0;
    
    totalp = s->branchp + s->shootp + s->stemp + s->structsurfp;
    s->inorglabp += totalp / 2.0;

    /* re-establish everything with C/N ~ 25 and C/P ~ 2500 */
    if (c->alloc_model == GRASSES) {
        s->branch = 0.0;
        s->branchn = 0.0;
        s->branchp = 0.0;
        s->sapwood = 0.0;
        s->stem = 0.0;
        s->stemn = 0.0;
        s->stemp = 0.0;
        s->stemnimm = 0.0;
        s->stemnmob = 0.0;
        s->stempimm = 0.0;
        s->stempmob = 0.0;
    } else {
        s->branch = 0.001;
        s->branchn = 0.00004;
        s->branchp = 0.0000004;
        s->sapwood = 0.001;
        s->stem = 0.001;
        s->stemn = 0.00004;
        s->stemnimm = 0.00004;
        s->stemnmob = 0.0;
        s->stemp = 0.0000004;
        s->stempimm = 0.0000004;
        s->stempmob = 0.0;
    }

    s->age = 0.0;
    s->metabsurf = 0.0;
    s->metabsurfn = 0.0;
    s->metabsurfp = 0.0;
    s->prev_sma = 1.0;
    s->root = 0.001;
    s->rootn = 0.00004;
    s->rootp = 0.0000004;
    s->shoot = 0.01;

    s->lai = p->sla * M2_AS_HA / KG_AS_TONNES / p->cfracts * s->shoot;
    s->shootn = 0.004;
    s->shootp = 0.00004;
    s->structsurf = 0.001;
    s->structsurfn = 0.00004;
    s->structsurfp = 0.0000004;

    /* reset litter flows */
    f->deadroots = 0.0;
    f->deadstems = 0.0;
    f->deadbranch = 0.0;
    f->deadsapwood = 0.0;
    f->deadleafn = 0.0;
    f->deadrootn = 0.0;
    f->deadbranchn = 0.0;
    f->deadstemn = 0.0;
    f->deadleafp = 0.0;
    f->deadrootp = 0.0;
    f->deadbranchp = 0.0;
    f->deadstemp = 0.0;

    /* update N:C, P:C of plant pools */
    if (float_eq(s->shoot, 0.0)) {
        s->shootnc = 0.0;
        s->shootpc = 0.0;
    } else {
        s->shootnc = s->shootn / s->shoot;
        s->shootpc = s->shootp / s->shoot;
    }

    if (c->ncycle == FALSE)
        s->shootnc = p->prescribed_leaf_NC;
    
    if (c->pcycle == FALSE)
        s->shootpc = p->prescribed_leaf_PC;

    if (float_eq(s->root, 0.0)) {
        s->rootnc = 0.0;
        s->rootpc = 0.0;
    } else {
        s->rootnc = MAX(0.0, s->rootn / s->root);
        s->rootpc = MAX(0.0, s->rootp / s->root);
    }
    
    return;
}

void hurricane(fluxes *f, params *p, state *s) {
    /* Specifically for the florida simulations - reduce LAI by 40%  */

    double orig_shoot_c, lost_c, lost_n, lost_p;
    double nc_leaf_litter, pc_leaf_litter;
    double lnleaf, fmleaf;
    ;
    /* Reduce LAI by 40%  */
    s->lai -= s->lai * 0.4;

    /* adjust C in the foliage */
    orig_shoot_c = s->shoot;
    s->shoot = s->lai / (p->sla * M2_AS_HA / KG_AS_TONNES / p->cfracts);
    lost_c = orig_shoot_c - s->shoot;
    lost_n = s->shootnc * lost_c;
    lost_p = s->shootpc * lost_c;
    s->shootn -= lost_n;
    s->shootp -= lost_p;

    /* Drop straight to floor, no retranslocation */

    /* C -> structural */
    if (float_eq(lost_c, 0.0)) {
        nc_leaf_litter = 0.0;
        pc_leaf_litter = 0.0;
    } else {
        nc_leaf_litter = lost_n / lost_c;
        pc_leaf_litter = lost_p / lost_c;
    }

    if (float_eq(nc_leaf_litter, 0.0)) {
        /* catch divide by zero if we have no leaves  */
        lnleaf = 0.0;
    } else {
        lnleaf = p->ligshoot / p->cfracts / nc_leaf_litter;
    }

    fmleaf = MAX(0.0, 0.85 - (0.018 * lnleaf));
    f->surf_struct_litter += lost_c * (1.0 - fmleaf);

    /* C -> metabolic */
    f->surf_metab_litter += lost_c * fmleaf;

    /* N -> structural */
    if (float_eq(f->surf_struct_litter, 0.0)) {
        f->n_surf_struct_litter += 0.0;
    } else {
        f->n_surf_struct_litter += (lost_n * f->surf_struct_litter *
                                    p->structrat / f->surf_struct_litter);
    }
    
    /* P -> structural */
    if (float_eq(f->surf_struct_litter, 0.0)) {
        f->p_surf_struct_litter += 0.0;
    } else {
        f->p_surf_struct_litter += (lost_p * f->surf_struct_litter *
                                    p->structratp / f->surf_struct_litter);
    }

    /* N -> metabolic pools */
    f->n_surf_metab_litter += lost_n - f->n_surf_struct_litter;
    
    /* P -> metabolic pools */
    f->p_surf_metab_litter += lost_p - f->p_surf_struct_litter;
    
    /* s->structsurf += lost_c; */
    /* s->structsurfn += lost_n; */

    
    return;
}
