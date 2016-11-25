#include "phenology.h"

void phenology(control *c, fluxes *f, met_arrays *ma, params *p, state *s) {
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
    int leaf_on_found, leaf_off_found;
    int project_day = c->day_idx;
    double grass_temp_threshold, tmax_ann, Tmin_avg, ppt_sum_crit;
    double gdd_thresh=-999.9;


    /*
        Krinner et al. 2005, page 26, alternatively Foley et al. 1996 suggests
        the same value = 100 for both pathways
    */
    if (c->alloc_model == GRASSES) {
        if (c->ps_pathway == C3)
            gdd_thresh = 185.;
        else if (c->ps_pathway == C4)
            gdd_thresh = 400.;
        calc_ini_grass_pheno_stuff(c, ma, project_day, &grass_temp_threshold,
                                   &tmax_ann, &Tmin_avg, &ppt_sum_crit);
    } else {
        gdd_thresh = gdd_chill_thresh(pa, pb, pc, p->previous_ncd);
    }


    calculate_leafon_off(c, ma, p, s, grass_temp_threshold, tmax_ann,
                         Tmin_avg, ppt_sum_crit, project_day,
                         &leaf_on, &leaf_off, &leaf_on_found,
                         &leaf_off_found, gdd_thresh);

    /*
        No leaf drop found, try a warmer temperature i.e. 5 instead of 0,
        if this doesn't work there really is an issue (or there is no leaf
        drop and we have an evergreen grass...
    */
    if (leaf_off_found == FALSE) {
        grass_temp_threshold = 5.0;
        calculate_leafon_off(c, ma, p, s, grass_temp_threshold, tmax_ann,
                             Tmin_avg, ppt_sum_crit, project_day,
                             &leaf_on, &leaf_off, &leaf_on_found,
                             &leaf_off_found, gdd_thresh);
    }

    /*
        if widening the temperature threshold didn't produce a suitable leaf
        drop date we will follow biome-bgc and assume the leaves fall on the
        last day
    */
    if (leaf_off_found == FALSE) {
        leaf_off = 364;
    }



    if (leaf_on_found == FALSE) {
        fprintf(stderr, "Problem in phenology leaf *ON* not found\n");
        exit(EXIT_FAILURE);
    }


    /*
        Length of time taken for new growth from storage to be allocated.
        This is either some site-specific calibration or the midpoint of the
        length of the growing season. The litterfall takes place over an
        identical period. Dividing by a larger number would increase the
        rate the C&N is allocated.
    */
    p->growing_seas_len = leaf_off - leaf_on;
    if (p->store_transfer_len < -900)
        len_groloss = (int)floor((float)p->growing_seas_len / 2.0);
    else
        len_groloss = p->store_transfer_len;

    calculate_days_left_in_growing_season(c, s, leaf_on, leaf_off, len_groloss);
    calculate_growing_season_fluxes(f, s, len_groloss);

    /*printf("%d %d\n", leaf_on, leaf_off); */

    return;
}

void calculate_leafon_off(control *c, met_arrays *ma, params *p, state *s,
                          double grass_temp_threshold, double tmax_ann,
                          double Tmin_avg, double ppt_sum_crit,
                          int project_day, int *leaf_on, int *leaf_off,
                          int *leaf_on_found, int *leaf_off_found,
                          double gdd_thresh) {

    double ppt_sum_next, ppt_sum, ppt_sum_prev, Tmean, Tsoil, Tsoil_next_3days,
           Tday;
    /*double Tmax; */
    double accumulated_ncd = 0.0;
    double accum_gdd = 0.0;
    /*double Tmin_boxcar; */
    int    drop_leaves = FALSE;
    int    d, dd, st, en, nov_doy;


    *leaf_on_found = FALSE;
    *leaf_off_found = FALSE;

    if (c->num_days == 366)
        nov_doy = 306;
    else
        nov_doy = 305;

    ppt_sum = 0.0;
    for (d = 1; d < c->num_days+1; d++) {
        Tmean = ma->tair[project_day];
        Tday = ma->tday[project_day];
        Tsoil = ma->tsoil[project_day];
        /*Tmax = ma->tmax[project_day];*/
        ppt_sum += ma->rain[project_day];

        /* Calculate ppt total from the next 7 days */
        if (d < 358) {
            st = project_day + 1;
            en = project_day + 8;
            ppt_sum_next = 0.0;
            for (dd = st; dd < en; dd++) {
                ppt_sum_next += ma->rain[dd];
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
                ppt_sum_prev += ma->rain[dd];
            }
        }

        if (d < 362) {
            Tsoil_next_3days = ((ma->tsoil[project_day] +
                                 ma->tsoil[project_day+1] +
                                 ma->tsoil[project_day+2]) / 3.0);

            /*Tmin_boxcar = ((ma->tmin[project_day-1] +
                            ma->tmin[project_day] +
                            ma->tmin[project_day+1]) / 3.0);*/

        } else {
            /* i.e. end of year, didn't find this so have no effect */
            Tsoil_next_3days = 999.9;
        }

        /* Sum the daily mean air temperature above 5degC starting on Jan 1 */
        accum_gdd += calc_gdd(Tmean);

        /*
        ** Calculate leaf on
        */
        if (c->alloc_model == GRASSES) {
            if (*leaf_on_found == FALSE &&
                accum_gdd >= gdd_thresh &&
                ppt_sum >= ppt_sum_crit) {

                *leaf_on = d;
                *leaf_on_found = TRUE;
            }
        } else {
            if (*leaf_on_found == FALSE && accum_gdd >= gdd_thresh) {
                  *leaf_on = d;
                  *leaf_on_found = TRUE;
            }
        }

        /*
        ** Calculate leaf off
        */
        if (c->alloc_model == GRASSES) {
            if (*leaf_off_found == FALSE) {

                /*
                    Leaf drop constraint is based on Foley et al. 1996

                    The 243 is just a safe guard to make sure we avoid
                    predicting offset in late spring (White et al. 1997).

                    This Tmean is the mean daytime temp, but I wonder if it
                    should be the full 24 daytime temp mean?
                */

                /*if (d >= 243) {
                    printf("%d %f %f\n", d, Tmean, grass_temp_threshold);
                }*/
                if (d >= 243 && Tday <= grass_temp_threshold) {
                    *leaf_off_found = TRUE;
                    *leaf_off = d;
                }


                /*
                    Leaf drop constraint is based on White et al. 1997

                     - test for hot && dry conditions.

                if (ppt_sum_prev < 11.4 &&
                    ppt_sum_next < 9.7 &&
                    Tmax > tmax_ann &&
                    d > 243) {

                    *leaf_off_found = TRUE;
                    *leaf_off = d;

                    - test for cold offset condition
                } else if (d > 243 && Tmin_boxcar <= Tmin_avg) {

                    *leaf_off_found = TRUE;
                    *leaf_off = d;

                }
                */

            }
        } else {
            if (*leaf_off_found == FALSE && accum_gdd >= gdd_thresh) {
                /*
                    I am prescribing that no leaves can fall off before doy=180
                    Had issue with KSCO simulations where the photoperiod was
                    less than the threshold very soon after leaf out.
                */
                if (d > 182) {
                    drop_leaves = leaf_drop(s->day_length[d-1], Tsoil,
                                            Tsoil_next_3days);
                    if (drop_leaves) {
                        *leaf_off_found = TRUE;
                        *leaf_off = d;
                    }
                }
            }
        }
        /* Calculated NCD from fixed date following Murray et al 1989. */
        if (d+1 >= nov_doy) {
            accumulated_ncd += calc_ncd(Tmean);
        }
        project_day++;
    }


    /* updated stored param, note this will be written out if the user
       dumps the current state, which makes sense as we may want pass the
       stat between spinup and a simulation */
    p->previous_ncd = accumulated_ncd;



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

void calc_ini_grass_pheno_stuff(control *c, met_arrays *ma, int project_day,
                                double *grass_temp_threshold,
                                double *tmax_ann, double *Tmin_avg,
                                double *ppt_sum_crit) {
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
    *Tmin_avg = 0.0;
    double tmin_ann = 70.0;
    double tavg_ann = 0.0;
    double ppt_sum = 0.0;
    double tair, tmin, tmax, Trange;

    for (d = 0; d < c->num_days; d++) {
        tair = ma->tair[project_day];
        tmax = ma->tmax[project_day];
        tmin = ma->tmin[project_day];
        ppt_sum += ma->rain[project_day];
        *Tmin_avg += ma->tmin[project_day];

        if (tmax > *tmax_ann)
           *tmax_ann = tmax;

        if (tmin < tmin_ann)
           tmin_ann = tmin;

        tavg_ann += tair;
        project_day += 1;
    }
    *Tmin_avg /= (float)c->num_days;

    /* reset date index */
    project_day = project_day_save;

    Trange = *tmax_ann - tmin_ann;
    tavg_ann /= c->num_days;

    /*
        Cool or warm grassland Definitions are from Botta, Table 1, pg 712.
        But grass temp thresholds are from Foley et al.
    */

    /* cool */
    if (Trange > 20.0 || tmin_ann < 5.0)
        *grass_temp_threshold = 0.0;

    /* warm */
    else if (Trange <= 20.0 || tmin_ann >= 5.0)
        *grass_temp_threshold = 5.0;
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

void calculate_growing_season_fluxes(fluxes *f, state *s, int len_groloss) {

    double denominator = (double)(len_groloss * len_groloss);

    /* C allocation rates across growing season */
    f->lrate = 2.0 * s->c_to_alloc_shoot / denominator;
    f->wrate = 2.0 * s->c_to_alloc_stem / denominator;
    f->brate = 2.0 * s->c_to_alloc_branch / denominator;
    f->crate = 2.0 * s->c_to_alloc_croot / denominator;

    /* N allocation rates across growing season */
    f->lnrate = 2.0 * s->n_to_alloc_shoot / denominator;
    f->bnrate = 2.0 * s->n_to_alloc_branch / denominator;
    f->wnimrate = 2.0 * s->n_to_alloc_stemimm / denominator;
    f->wnmobrate = 2.0 * s->n_to_alloc_stemmob / denominator;
    f->cnrate = 2.0 * s->n_to_alloc_croot / denominator;
    
    /* P allocation rates across growing season */
    f->lprate = 2.0 * s->p_to_alloc_shoot / denominator;
    f->bprate = 2.0 * s->p_to_alloc_branch / denominator;
    f->wpimrate = 2.0 * s->p_to_alloc_stemimm / denominator;
    f->wpmobrate = 2.0 * s->p_to_alloc_stemmob / denominator;
    f->cprate = 2.0 * s->p_to_alloc_croot / denominator;

    return;
}
