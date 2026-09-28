#include "disturbance.h"


void figure_out_years_with_disturbances(control *c, met_arrays *ma, params *p,
                                        int **yrs, int *cnt) {
    /*
        Years in which a fire happens (on disturbance_doy). Either a single
        prescribed year (burn_specific_yr) or every return_interval years,
        starting return_interval years after the first year of the forcing.
    */
    int first_year, last_year, year;

    /* the met arrays are in timestep order, so the ends give the years */
    if (c->sub_daily) {
        first_year = (int)ma->year[0];
        last_year = (int)ma->year[c->total_num_days * c->num_hlf_hrs - 1];
    } else {
        first_year = (int)ma->year[0];
        last_year = (int)ma->year[c->total_num_days - 1];
    }

    *cnt = 0;
    if (p->burn_specific_yr > -900) {
        (*yrs)[0] = p->burn_specific_yr;
        *cnt = 1;
    } else {
        for (year = first_year + time_till_next_disturbance(p);
             year <= last_year;
             year += time_till_next_disturbance(p)) {

            if (*cnt > 0) {
                *yrs = (int *)realloc(*yrs, (*cnt + 1) * sizeof(int));
                if (*yrs == NULL) {
                    fprintf(stderr,"Error resizing years array\n");
                    exit(EXIT_FAILURE);
                }
            }
            (*yrs)[*cnt] = year;
            *cnt += 1;
        }
    }

    return;
}

int time_till_next_disturbance(params *p) {
    /* calculate the number of years until a disturbance event occurs
    assuming a return interval of X years. Deterministic for now, a random
    interval would be -log(1 - U) * return_interval (Knuth 3.4.1). */

    return (MAX(1, p->return_interval));
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
    * 50 percent of N returned to inorgn pool"
    * P isn't volatilised: all the burnt P is returned to the labile pool
      (the sorption re-equilibrates the next day)
    * Coarse roots are not damaged by fire!

    vaguely following ...
    http://treephys.oxfordjournals.org/content/24/7/765.full.pdf
    */
    double totaln;

    totaln = s->branchn + s->shootn + s->stemn + s->structsurfn;
    s->inorgn += totaln / 2.0;
    if (c->pcycle) {
        s->inorglabp += s->branchp + s->shootp + s->stemp + s->structsurfp;
    }

    /* re-establish everything with C/N ~ 25.  */
    if (c->alloc_model == GRASSES) {
        s->branch = 0.0;
        s->branchn = 0.0;
        s->sapwood = 0.0;
        s->stem = 0.0;
        s->stemn = 0.0;
        s->stemnimm = 0.0;
        s->stemnmob = 0.0;
    } else {
        s->branch = 0.001;
        s->branchn = 0.00004;
        s->sapwood = 0.001;
        s->stem = 0.001;
        s->stemn = 0.00004;
        s->stemnimm = 0.00004;
        s->stemnmob = 0.0;
    }

    s->age = 0.0;
    s->metabsurf = 0.0;
    s->metabsurfn = 0.0;
    s->prev_sma = 1.0;
    s->root = 0.001;
    s->rootn = 0.00004;
    s->shoot = 0.01;

    s->lai = p->sla * M2_AS_HA / KG_AS_TONNES / p->cfracts * s->shoot;
    s->shootn = 0.004;
    s->structsurf = 0.001;
    s->structsurfn = 0.00004;

    /* reset litter flows */
    f->deadroots = 0.0;
    f->deadstems = 0.0;
    f->deadbranch = 0.0;
    f->deadsapwood = 0.0;
    f->deadleafn = 0.0;
    f->deadrootn = 0.0;
    f->deadbranchn = 0.0;
    f->deadstemn = 0.0;

    /* P pools re-established at a N:P of 10 */
    if (c->pcycle) {
        s->branchp = s->branchn / 10.0;
        s->stempimm = s->stemnimm / 10.0;
        s->stempmob = s->stemnmob / 10.0;
        s->stemp = s->stempimm + s->stempmob;
        s->metabsurfp = 0.0;
        s->rootp = s->rootn / 10.0;
        s->shootp = s->shootn / 10.0;
        s->structsurfp = s->structsurfn / 10.0;
        f->deadleafp = 0.0;
        f->deadrootp = 0.0;
        f->deadbranchp = 0.0;
        f->deadstemp = 0.0;
        s->shootpc = s->shootp / s->shoot;
        s->rootpc = s->rootp / s->root;
    }

    /* update N:C of plant pools */
    if (float_eq(s->shoot, 0.0))
        s->shootnc = 0.0;
    else
        s->shootnc = s->shootn / s->shoot;

    if (c->ncycle == FALSE)
        s->shootnc = p->prescribed_leaf_NC;

    if (float_eq(s->root, 0.0))
        s->rootnc = 0.0;
    else
        s->rootnc = MAX(0.0, s->rootn / s->root);

    return;
}

void hurricane(fluxes *f, params *p, state *s) {
    /*
        Specifically for the florida simulations - reduce LAI by 40%

        Called after calculate_litterfall, so the lost foliage is added to
        the day's leaf litter; update_plant_state then removes the C from the
        shoot, carbon_allocation reduces the LAI accordingly and the litter
        is partitioned into the soil pools as usual. Previously the lost C&N
        were added to the partitioned litter fluxes, which are overwritten
        later in the day, i.e. they vanished. No retranslocation.
    */
    double frac_lost = 0.4, lost_c, lost_n, lost_p;

    lost_c = s->shoot * frac_lost;
    lost_n = s->shootn * frac_lost;
    lost_p = s->shootp * frac_lost;

    f->deadleaves += lost_c;
    f->deadleafn += lost_n;
    f->deadleafp += lost_p;

    /* shoot N & P aren't updated from the dead leaf N & P, so remove them
       here */
    s->shootn -= lost_n;
    s->shootp -= lost_p;

    return;
}
