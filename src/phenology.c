#include "phenology.h"

void phenology(control *c, fluxes *f, met *m, params *p, state *s,
               double *daylen, int project_day) {
    /*
    There are two phenology schemes currently implemented, one which should
    generally be applicable for deciduous broadleaf forests && one for
    short-lived grasses.

    The tree scheme is based on Botta et al. who calibrated their model against
    EO data. The grasses scheme is really a combination of approaches, I've
    tried to keep it is as simple as possible...user beware :)

    There are some key issues:
     - leaf drop for trees won't work where the daylength isn't applicable, i.e
       the tropics would be my guess.
     - the grass phenology drop takes no account of available water. It of
       course should, but this would need a complete re-write of the logic.
       Perhaps someone else can do that.

    One potential thing to look at if you ever get bored is:

     - Look at Caldararu et al 2014, Phenology as a strategy for carbon
       optimality: a global model. Seems interesting at a quick glance.

    The distribution of C&N is pre-calculated here using a ramping function
    based on leaf out/off dates.

    Finally, no account has been taken for the southern hemisphere! This won't
    work there.

    References:
    -----------
    * Botta, A. et al. (2000) GCB, 6, 709-725.
    * Foley, J. et al. (1996) GBC, 10, 603-628.
    * Krinner, G. et al. (2005) GBC, 19, GB1015
    * White, M. A. et al. (1997) GBC, 11, 217-234.
    */

    /* (days) Leaf flush params following Botta. */
    double pa = -68.0;

    /* (days) Leaf flush params following Botta. */
    double pb = 638.0;

    /* (1/days) Leaf flush params following Botta. */
    double pc = -0.01;

    int leaf_on = 0, leaf_off = 0, len_groloss = 0.0;

    calculate_leafon_off(c, m, p, daylen, project_day, pa, pb, pc, &leaf_on,
                         &leaf_off, &len_groloss);

    calculate_days_left_in_growing_season(c, s, leaf_on, leaf_off, len_groloss);
    calculate_growing_season_fluxes(f, s, len_groloss);

    return;
}

void calculate_leafon_off(control *c, met *m, params *p, double *daylen,
                         int project_day, double pa, double pb, double pc,
                         int *leaf_on, int *leaf_off, int *len_groloss) {

    double grass_temp_threshold, tmax_ann, ppt_sum_crit, ppt_sum_next, ppt_sum,
           ppt_sum_prev, Tmean, Tsoil, Tsoil_next_3days, Tair_next_3days,
           gdd_thresh;
    double accumulated_ncd = 0.0;
    double accum_gdd = 0.0;
    int    leaf_on_found = FALSE;
    int    leaf_off_found = FALSE;
    int    drop_leaves = FALSE;
    int    d, dd, st, en, nov_doy;

    /*
        Krinner et al. 2005, page 26, alternatively Foley et al. 1996 suggests
        the same value = 100 for both pathways
    */
    if (c->alloc_model == GRASSES) {
        if (c->ps_pathway == C3)
            gdd_thresh = 185.;
        else if (c->ps_pathway == C4)
            gdd_thresh = 400.;
        calc_ini_grass_pheno_stuff(c, m, project_day, &grass_temp_threshold,
                                   &tmax_ann, &ppt_sum_crit);
    } else {
        gdd_thresh = gdd_chill_thresh(pa, pb, pc, p->previous_ncd);
    }

    if (c->num_days == 366)
        nov_doy = 306;
    else
        nov_doy = 305;

    ppt_sum = 0.0;
    for (d = 1; d < c->num_days+1; d++) {
        Tmean = m->tair[project_day];
        Tsoil = m->tsoil[project_day];
        ppt_sum += m->rain[project_day];

        /* Calculate ppt total from the next 7 days */
        if (d < 358) {
            st = project_day + 1;
            en = project_day + 8;
            ppt_sum_next = 0.0;
            for (dd = st; dd < en; dd++) {
                ppt_sum_next += m->rain[dd];
            }
        } else {
            /* i.e. end of year, didn't find this so have no effect */
            ppt_sum_next = 0.0;
        }

        /* Calculate ppt total from the previous 30 days */
        if (project_day < 30) {
            ppt_sum_prev = 0.0;
        } else {
            st = project_day - 30;
            en = project_day;
            ppt_sum_prev = 0.0;
            for (dd = st; dd < en; dd++) {
                ppt_sum_prev += m->rain[dd];
            }
        }

        if (d < 362) {
            Tsoil_next_3days = ((m->tsoil[project_day] +
                                 m->tsoil[project_day+1] +
                                 m->tsoil[project_day+2]) / 3.0);

            Tair_next_3days = ((m->tair[project_day] +
                                m->tair[project_day+1] +
                                m->tair[project_day+2]) / 3.0);
        } else {
            /* i.e. end of year, didn't find this so have no effect */
            Tsoil_next_3days = 999.9;
            Tair_next_3days = 999.9;
        }

        /* Sum the daily mean air temperature above 5degC starting on Jan 1 */
        accum_gdd += calc_gdd(Tmean);

        /*
        ** Calculate leaf on
        */
        if (c->alloc_model == GRASSES) {
            if (leaf_on_found == FALSE &&
                accum_gdd >= gdd_thresh &&
                ppt_sum >= ppt_sum_crit) {

                *leaf_on = d;
                leaf_on_found = TRUE;
            }
        } else {
            if (leaf_on_found == FALSE && accum_gdd >= gdd_thresh) {
                  *leaf_on = d;
                  leaf_on_found = TRUE;
            }
        }

        /*
        ** Calculate leaf off
        */
        if (c->alloc_model == GRASSES) {
            if (leaf_off_found == FALSE) {
                /* test for hot && dry conditions, based on white et al. 1997 */
                if (ppt_sum_prev < 11.4 &&
                    ppt_sum_next < 9.7 &&
                    Tmean > tmax_ann) {

                    leaf_off_found = TRUE;
                    *leaf_off = d;
                } else if (d > 182 && Tair_next_3days < grass_temp_threshold) {
                    /*
                        test for cold offset condition
                        Leaf drop constraint is based on Foley et al. 1996 as
                        we dont have access to the Tmin for the constraint from
                        White et al. This is likely more straightforward anyway
                    */
                    leaf_off_found = TRUE;
                    *leaf_off = d;
                }
            }
        } else {
            if (leaf_off_found == FALSE && accum_gdd >= gdd_thresh) {
                /*
                    I am prescribing that no leaves can fall off before doy=180
                    Had issue with KSCO simulations where the photoperiod was
                    less than the threshold very soon after leaf out.
                */
                if (d > 182) {
                    drop_leaves = leaf_drop(daylen[d-1], Tsoil, Tsoil_next_3days);
                    if (drop_leaves) {
                        leaf_off_found = TRUE;
                        *leaf_off = d;
                    }
                }
            }
        }
        /* Calculated NCD from fixed date following Murray et al 1989. */
        if (d+1 >= nov_doy)
            accumulated_ncd += calc_ncd(Tmean);
        project_day++;
    }

    /* updated stored param, note this will be written out if the user
       dumps the current state, which makes sense as we may want pass the
       stat between spinup and a simulation */
    p->previous_ncd = accumulated_ncd;

    /*
        Length of time taken for new growth from storage to be allocated.
        This is either some site-specific calibration or the midpoint of the
        length of the growing season. The litterfall takes place over an
        identical period. Dividing by a larger number would increase the
        rate the C&N is allocated.
    */
    p->growing_seas_len = *leaf_off - *leaf_on;
    if (p->store_transfer_len < -900)
        *len_groloss = (int)floor((float)p->growing_seas_len / 2.0);
    else
        *len_groloss = p->store_transfer_len;

    if (leaf_on_found == FALSE || leaf_off_found == FALSE) {
        fprintf(stderr, "Problem in phenology leaf on/off not found\n");
        exit(EXIT_FAILURE);
    }
    return;
}


double calc_gdd(double Tavg) {
    /*
        calculate the number of growing degree days, hypothesis is that
        leaves appear after a threshold has been passed.
    */
    double Tbase = 5.0; /* degC */

    return (MAX(0.0, Tavg - Tbase));
}

double gdd_chill_thresh(double pa, double pb, double pc, double ncd) {
    /*
        Leaf out has a chilling requirement, num chill days reduces the GDD
        requirement
    */
    return (pa + pb * exp(pc * ncd));
}

double calc_ncd(double Tmean) {
    /*
        Calculate the number of chilling days from fixed dates (1 Nov),
        following Murray et al. 1989, same as Botta does. """
    */
    if (Tmean < 5.0)
        return (1.0);
    else
        return (0.0);
}

double leaf_drop(double daylen, double Tsoil, double Tsoil_next_3days) {
    /*
        Thresholds to drop leaves come from White et al.
        Note 655 minutes = 10.916 hrs.

        - Dependance on daylength means that this is only valid outside of the
          tropics.

        References:
        -----------
        White, M. A. et al. (1997) GBC, 11, 217-234.
    */
    if ((daylen <= 10.9166667 && Tsoil <= 11.15) || Tsoil_next_3days < 2.0)
        return (TRUE);
    else
        return (FALSE);
}

void calc_ini_grass_pheno_stuff(control *c, met *m, int project_day,
                                double *grass_temp_threshold,
                                double *tmax_ann, double *ppt_sum_crit) {
    /*
        Series of constraints based on temp && precip need to be
        pre-calculated for grasses to determine leaf on/off
    */

    /*
        Save this as we need to loop over the data once to pre-calculate
        everything
    */
    int project_day_save = project_day, d;
    *tmax_ann = 0.0;
    double tmin_ann = 70.0;
    double tavg_ann = 0.0;
    double ppt_sum = 0.0;
    double tair, tam, tpm, Trange;

    for (d = 0; d < c->num_days; d++) {
        tair = m->tair[project_day];
        tam = m->tam[project_day];
        tpm = m->tpm[project_day];
        ppt_sum += m->rain[project_day];

        if (tair > *tmax_ann)
           *tmax_ann = tair;

        if (tair < tmin_ann)
           tmin_ann = tair;

        tavg_ann += tair;
        project_day += 1;
    }

    /* reset date index */
    project_day = project_day_save;

    Trange = *tmax_ann - tmin_ann;
    tavg_ann /= c->num_days;

    /*
        Cool or warm grassland Definitions are from Botta, Table 1, pg 712.
        But thresholds are from Foley et al.
    */
    if (Trange > 20.0 || tmin_ann < 5.0)
        *grass_temp_threshold = 0.0; /* cool */
    else if (Trange <= 20.0 || tmin_ann >= 5.0)
        *grass_temp_threshold = 5.0;  /* warm */
    else {
        fprintf(stderr, "Problem grass thresholds\n");
        exit(EXIT_FAILURE);
    }

    /*
        92% of tmax_ann is the threshold used in grass offset below
        Note this has to be done below the range calcs as they use the tmax
    */
    *tmax_ann *= 0.92;

    /*
        Based on White et al. 1997 this is a threshold for grasses so they
        have enough accumulated rain. It is essentially a fudge for soil
        moisture availability && the 15% is somewhat arbitary
    */
    *ppt_sum_crit = ppt_sum * 0.15;

    return;
}

void calculate_days_left_in_growing_season(control *c, state *s,
                                            int leaf_on, int leaf_off,
                                            int len_groloss) {
    /* Calculate 2 arrays to signify the days left of growing period
    an days left before all the leaves fall off. In both cases these will
    be 2 lists, with 0.0 outside of the growing period and a series of
    numbers e.g. day 46 to 0 for growing_days. 0.5 is subtracted from the
    doy to get round the issue of approximating an integral with discrete
    time steps -> trapezoidal type solution
    */
    int doy;

    for (doy = 1; doy < c->num_days+1; doy++) {

        if (doy > leaf_off - len_groloss && doy <= leaf_off) {
            s->remaining_days[doy-1] = (doy - 0.5) - leaf_off + len_groloss;
        } else {
            s->remaining_days[doy-1] = 0.0;
        }

        if (doy > leaf_on && doy <= len_groloss+leaf_on) {
            s->growing_days[doy-1] = len_groloss + leaf_on - (doy - 0.5);
        } else {
            s->growing_days[doy-1] = 0.0;
        }

        if (doy > leaf_on && doy < leaf_off) {
            s->leaf_out_days[doy-1] = 1.0;
        } else {
            s->leaf_out_days[doy-1] = 0.0;
        }
    }

    return;
}

void calculate_growing_season_fluxes(fluxes *f, state *s, double len_groloss) {

    /* C allocation rates across growing season */
    f->lrate = 2.0 * s->c_to_alloc_shoot / (len_groloss * len_groloss);
    f->wrate = 2.0 * s->c_to_alloc_stem / (len_groloss * len_groloss);
    f->brate = 2.0 * s->c_to_alloc_branch / (len_groloss * len_groloss);
    f->crate = 2.0 * s->c_to_alloc_croot / (len_groloss * len_groloss);

    /* N allocation rates across growing season */
    f->lnrate = 2.0 * s->n_to_alloc_shoot / (len_groloss * len_groloss);
    f->bnrate = 2.0 * s->n_to_alloc_branch / (len_groloss * len_groloss);
    f->wnimrate = 2.0 * s->n_to_alloc_stemimm / (len_groloss * len_groloss);
    f->wnmobrate = 2.0 * s->n_to_alloc_stemmob / (len_groloss * len_groloss);
    f->cnrate = 2.0 * s->n_to_alloc_croot / (len_groloss * len_groloss);

    return;
}
