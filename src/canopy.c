/* ============================================================================
* Calculates all within canopy C & water fluxes (live in water balance).
*
*
* NOTES:
*   - Should restructure the code so that MATE is called from within the canopy
*     space, rather than via plant growth
*
*   Future improvements:
*    - Add a two-stream approximation.
*    - Add a clumping term to the extinction coefficients for apar calcs
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   09.02.2016
*
* =========================================================================== */
#include "canopy.h"

void canopy(canopy_wk *cw, control *c, fluxes *f, met_arrays *ma, met *m,
            nrutil *nr, params *p, state *s) {
    /*
        Canopy module consists of two parts:
        (1) a radiation sub-model to calculate apar of sunlit/shaded leaves
            - this is all handled in radiation.c
        (2) a coupled model of stomatal conductance, photosynthesis and
            the leaf energy balance to solve the leaf temperature and partition
            absorbed net radiation between sensible and latent heat.
        - The canopy is represented by a single layer with two big leaves
          (sunlit & shaded).

        - The logic broadly follows MAESTRA code, with some restructuring.

        References
        ----------
        * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.
        * Dai et al. (2004) Journal of Climate, 17, 2281-2299.
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
    */
    int    hod, iter = 0, itermax = 100, dummy=0, sunlight_hrs;
    int    debug = TRUE, stressed = FALSE;
    double doy, year, dummy2=0.0, previous_sw, current_sw, gsv;
    double previous_cs, current_cs;

    /* loop through the day */
    zero_carbon_day_fluxes(f);
    zero_water_day_fluxes(f);
    previous_sw = s->pawater_topsoil + s->pawater_root;
    previous_cs = s->canopy_store;
    sunlight_hrs = 0;
    doy = ma->doy[c->hour_idx];
    year = ma->year[c->hour_idx];

    for (hod = 0; hod < c->num_hlf_hrs; hod++) {
        unpack_met_data(c, f, ma, m, hod, dummy2);

        /* calculates diffuse frac from half-hourly incident radiation */
        calculate_solar_geometry(cw, p, doy, hod);
        get_diffuse_frac(cw, doy, m->sw_rad);
        
        unpack_solar_geometry(cw, c);

        /* Is the sun up? */
        if (cw->elevation > 0.0 && m->par > 20.0) {
            calculate_absorbed_radiation(cw, p, s, m->par);
            calculate_top_of_canopy_leafn(cw, p, s);
            calculate_top_of_canopy_leafp(cw, p, s);
            calc_leaf_to_canopy_scalar(cw, p);
            
            /* sunlit / shaded loop */
            for (cw->ileaf = 0; cw->ileaf < NUM_LEAVES; cw->ileaf++) {

                /* initialise values of Tleaf, Cs, dleaf at the leaf surface */
                initialise_leaf_surface(cw, m);

                /* Leaf temperature loop */
                while (TRUE) {

                    if (c->ps_pathway == C3) {
                        photosynthesis_C3(c, cw, m, p, s);
                    } else {
                        /* Nothing implemented */
                        fprintf(stderr, "C4 photosynthesis not implemented\n");
                        exit(EXIT_FAILURE);
                    }

                    if (cw->an_leaf[cw->ileaf] > 1E-04) {

                        if (c->water_balance == HYDRAULICS) {
                            // Ensure transpiration does not exceed Emax, if it
                            // does we recalculate gs and An
                            stressed = calculate_emax(c, cw, f, m, p, s);
                        }

                        if (stressed == FALSE) {
                            /* Calculate new Cs, dleaf, Tleaf */
                            solve_leaf_energy_balance(c, cw, f, m, p, s);
                        }
                    } else {
                        break;
                    }

                    if (iter >= itermax) {
                        fprintf(stderr, "No convergence in canopy loop:\n");
                        exit(EXIT_FAILURE);
                    } else if (fabs(cw->tleaf[cw->ileaf] - cw->tleaf_new) < 0.02) {
                        break;
                    }

                    /* Update temperature & do another iteration */
                    cw->tleaf[cw->ileaf] = cw->tleaf_new;
                    iter++;
                } /* end of leaf temperature loop */
            } /* end of sunlit/shaded leaf loop */
        } else {
            zero_hourly_fluxes(cw);

            /* set tleaf to tair during the night */
            cw->tleaf[SUNLIT] = m->tair;
            cw->tleaf[SHADED] = m->tair;

            /*
            ** pre-dawn soil water potential (MPa), clearly one should link this
            ** the actual sun-rise :). Here 10 = 5 am, 10 is num_half_hr
            **/
            if (c->water_balance == HYDRAULICS && hod == 10) {
                s->predawn_swp = s->weighted_swp;
                /*_calc_soil_water_potential(c, p, s);*/

            }

        }
        scale_leaf_to_canopy(c, cw, s);
        if (c->water_balance == HYDRAULICS && hod == 24) {
            s->midday_lwp = cw->lwp_canopy;
        }
        sum_hourly_carbon_fluxes(cw, f, p);

        calculate_water_balance_sub_daily(c, f, m, nr, p, s, dummy,
                                          cw->trans_canopy, cw->omega_canopy,
                                          cw->rnet_canopy);

        if (c->print_options == SUBDAILY && c->spin_up == FALSE) {
            write_subdaily_outputs_ascii(c, cw, year, doy, hod);
        }
        c->hour_idx++;
        sunlight_hrs++;
    } /* end of hour loop */

    /* work out average omega for the day over sunlight hours */
    f->omega /= sunlight_hrs;

    if (c->water_stress) {
        // Calculate the soil moisture availability factors [0,1] in the
        // topsoil and the entire root zone
        calculate_soil_water_fac(c, p, s);

        //printf("%lf %.10lf\n", s->wtfac_root, s->saved_swp);
    } else {
        s->wtfac_topsoil = 1.0;
        s->wtfac_root = 1.0;
    }

    current_sw = s->pawater_topsoil + s->pawater_root;
    current_cs = s->canopy_store;
    check_water_balance(c, f, s, previous_sw, current_sw, previous_cs,
                        current_cs, year, doy);

    //if (debug) {
    //    current_sw = s->pawater_topsoil + s->pawater_root;
    //    check_water_balance(c, f, previous_sw, current_sw);
    //}

    //if (c->pdebug) {
    //    exit(1);
    //}
    return;
}

void solve_leaf_energy_balance(control *c, canopy_wk *cw, fluxes *f, met *m,
                               params *p, state *s) {
    /*
        Wrapper to solve conductances, transpiration and calculate a new
        leaf temperautre, vpd and Cs at the leaf surface.

        - The logic broadly follows MAESTRA code, with some restructuring.

        References
        ----------
        * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.

    */
    int    idx;
    double omega, transpiration, LE, Tdiff, gv, gbc, gh, sw_rad;

    idx = cw->ileaf;
    sw_rad = cw->apar_leaf[idx] * PAR_2_SW; /* W m-2 */
    cw->rnet_leaf[idx] = calc_leaf_net_rad(p, s, m->tair, m->vpd, sw_rad);
    penman_leaf_wrapper(m, p, s, cw->tleaf[idx], cw->rnet_leaf[idx],
                        cw->gsc_leaf[idx], &transpiration, &LE, &gbc, &gh, &gv,
                        &omega);

    /* store in structure */
    cw->trans_leaf[idx] = transpiration;
    cw->omega_leaf[idx] = omega;

    /*
     * calculate new Cs, dleaf & tleaf
     */
    Tdiff = (cw->rnet_leaf[idx] - LE) / (CP * MASS_AIR * gh);
    cw->tleaf_new = m->tair + Tdiff / 4.0;
    cw->Cs = m->Ca - cw->an_leaf[idx] / gbc;
    cw->dleaf = cw->trans_leaf[idx] * m->press / gv;

    return;
}

double calc_leaf_net_rad(params *p, state *s, double tair, double vpd,
                         double sw_rad) {

    double rnet, Tk, ea, emissivity_atm, net_lw_rad;
    /*
        extinction coefficient for diffuse radiation and black leaves
        (m2 ground m2 leaf)
    */
    double kd = 0.8;

    /* isothermal net LW radiaiton at top of canopy, assuming emissivity of
       the canopy is 1 */
    Tk = tair + DEG_TO_KELVIN;

    /* Isothermal net radiation (Leuning et al. 1995, Appendix) */
    ea = calc_sat_water_vapour_press(tair) - vpd;

    /* catch for AWAP diurnal stuff until I better connect VPD and Tair */
    if (ea < 0.0) {
        ea = 0.0;
    }

    /* apparent emissivity for a hemisphere radiating at air temp eqn D4 */
    emissivity_atm = 0.642 * pow((ea / Tk), (1.0 / 7.0));

    net_lw_rad = (1.0 - emissivity_atm) * SIGMA * pow(Tk, 4.0);
    rnet = p->leaf_abs * sw_rad - net_lw_rad * kd * exp(-kd * s->lai);

    return (rnet);
}

void zero_carbon_day_fluxes(fluxes *f) {

    f->gpp_gCm2 = 0.0;
    f->npp_gCm2 = 0.0;
    f->gpp = 0.0;
    f->npp = 0.0;
    f->auto_resp = 0.0;
    f->apar = 0.0;

    return;
}


void calculate_top_of_canopy_leafn(canopy_wk *cw, params *p, state *s) {

    /*
    Calculate the N at the top of the canopy (g N m-2), N0.

    References:
    -----------
    * Chen et al 93, Oecologia, 93,63-69.

    */
    double Ntot;

    /* leaf mass per area (g C m-2 leaf) */
    double LMA = 1.0 / p->sla * p->cfracts * KG_AS_G;
    
    if (s->lai > 0.0) {
        /* the total amount of nitrogen in the canopy */
        Ntot = s->shootnc * LMA * s->lai;

        /* top of canopy leaf N (gN m-2) */
        cw->N0 = Ntot * p->kn / (1.0 - exp(-p->kn * s->lai));
    } else {
        cw->N0 = 0.0;
    }

    return;
}

void calculate_top_of_canopy_leafp(canopy_wk *cw, params *p, state *s) {
  
  /*
   Calculate the P at the top of the canopy (g P m-2), P0, based on
   calculate_top_of_canopy_leafp relationship;
   
   */
  double Ptot;
  
  /* leaf mass per area (g C m-2 leaf) */
  double LMA = 1.0 / p->sla * p->cfracts * KG_AS_G;
  
  if (s->lai > 0.0) {
    /* the total amount of phosphorus in the canopy */
    Ptot = s->shootpc * LMA * s->lai;
    
    /* top of canopy leaf N (gN m-2) */
    cw->P0 = Ptot * p->kp / (1.0 - exp(-p->kp * s->lai));
  } else {
    cw->P0 = 0.0;
  }
  
  return;
}

void zero_hourly_fluxes(canopy_wk *cw) {

    int i;

    /* sunlit / shaded loop */
    for (i = 0; i < NUM_LEAVES; i++) {
        cw->an_leaf[i] = 0.0;
        cw->rd_leaf[i] = 0.0;
        cw->gsc_leaf[i] = 0.0;
        cw->trans_leaf[i] = 0.0;
        cw->rnet_leaf[i] = 0.0;
        cw->apar_leaf[i] = 0.0;
        cw->omega_leaf[i] = 0.0;
    }

    return;
}

void scale_leaf_to_canopy(control *c, canopy_wk *cw, state *s) {

    double beta;

    cw->an_canopy = cw->an_leaf[SUNLIT] + cw->an_leaf[SHADED];
    cw->rd_canopy = cw->rd_leaf[SUNLIT] + cw->rd_leaf[SHADED];
    cw->gsc_canopy = cw->gsc_leaf[SUNLIT] + cw->gsc_leaf[SHADED];
    cw->apar_canopy = cw->apar_leaf[SUNLIT] + cw->apar_leaf[SHADED];
    cw->trans_canopy = cw->trans_leaf[SUNLIT] + cw->trans_leaf[SHADED];
    cw->omega_canopy = (cw->omega_leaf[SUNLIT] + cw->omega_leaf[SHADED]) / 2.0;
    cw->rnet_canopy = cw->rnet_leaf[SUNLIT] + cw->rnet_leaf[SHADED];

    if (c->water_balance == HYDRAULICS) {
        cw->lwp_canopy = (cw->lwp_leaf[SUNLIT] + cw->lwp_leaf[SHADED]) / 2.0;

        beta = (cw->fwsoil_leaf[SUNLIT] + cw->fwsoil_leaf[SHADED]) / 2.0;
        s->wtfac_topsoil = beta;
        s->wtfac_root = beta;
    }


    return;
}

void sum_hourly_carbon_fluxes(canopy_wk *cw, fluxes *f, params *p) {

    /* umol m-2 s-1 -> gC m-2 30 min-1 */
    f->gpp_gCm2 += cw->an_canopy * UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HLFHR;
    f->npp_gCm2 = f->gpp_gCm2 * p->cue;
    f->gpp = f->gpp_gCm2 * GRAM_C_2_TONNES_HA;
    f->npp = f->npp_gCm2 * GRAM_C_2_TONNES_HA;
    f->auto_resp = f->gpp - f->npp;

    /* umol m-2 s-1 -> J m-2 s-1 -> MJ m-2 30 min-1 */
    f->apar += cw->apar_canopy * UMOL_2_JOL * J_TO_MJ * SEC_2_HLFHR;
    f->gs_mol_m2_sec += cw->gsc_canopy;

    return;
}

void initialise_leaf_surface(canopy_wk *cw, met *m) {
    /* initialise values of Tleaf, Cs, dleaf at the leaf surface */
    cw->tleaf[cw->ileaf] = m->tair;
    cw->dleaf = m->vpd;
    cw->Cs = m->Ca;
}

void calc_leaf_to_canopy_scalar(canopy_wk *cw, params *p) {
    /*
        Calculate scalar to transform leaf Vcmax and Jmax values to big leaf
        values. Following Wang & Leuning, as long as sunlit and shaded
        leaves are treated seperately, values of parameters in the coupled
        model for the two big leaves can be closely approximated by
        integrating values for individual leaves.

        - Inserting eqn C6 & C7 into B5

        per unit ground area

        Parameters:
        ----------
        canopy_wk : structure
            various canopy values: in this case the sunlit or shaded LAI &
            cos_zenith angle.
        scalar_sun : float
            scalar for sunlit leaves, values returned in unit ground area
            (returned)
        scalar_sha : float
            scalar for shaded leaves, values returned in unit ground area
            (returned)

        References:
        ----------
        * Wang and Leuning (1998) AFm, 91, 89-111; particularly the Appendix.
    */
    double kn, lai_sun, lai_sha;
    kn = p->kn;
    lai_sun = cw->lai_leaf[SUNLIT];
    lai_sha = cw->lai_leaf[SHADED];

    cw->cscalar[SUNLIT] = (1.0 - exp(-(cw->kb + kn) * lai_sun)) / (cw->kb + kn);
    cw->cscalar[SHADED] = (1.0 - exp(-kn * lai_sha)) / kn - cw->cscalar[SUNLIT];

    return;
}

void check_water_balance(control *c, fluxes *f, state *s, double previous_sw,
                         double current_sw, double previous_cs,
                         double current_cs, double year, int doy) {

    double delta_sw, delta_cs;

    delta_sw = current_sw - previous_sw;
    delta_cs = current_cs - previous_cs;
    f->day_wbal = f->day_ppt - (f->runoff + f->et + delta_cs + delta_sw);

    printf("%d,%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
           (int)year, doy,
           f->day_ppt, f->runoff, \
           f->et, delta_sw, delta_cs, previous_sw, current_sw,
           s->predawn_swp, s->midday_lwp, f->day_wbal);


    return;
}

int calculate_emax(control *c, canopy_wk *cw, fluxes *f, met *m, params *p,
                    state *s) {

    // Assumption that during the day transpiration cannot exceed a maximum
    // value, Emax. At this point we've reached a leaf water potential minimum.
    // Once this point is reached transpiration, gs and A are reclulated
    //
    // Reference:
    // * Duursma et al. 2008, Tree Physiology 28, 265–276
    //

    // plant component of the leaf-specific hydraulic conductance
    // (mmol m–2 s–1 MPa–1)
    double kp = 2.0;
    double ktot, emax_leaf, etest, gsv, frac;
    int    idx = cw->ileaf;
    int    stressed = FALSE;

    // Need to work out the direct/diffuse weighting to adjust the
    // stressed gs/transpiration calculation.
    //if (cw->ileaf == 0) {
    //    frac = 1.0 - cw->diffuse_frac;
    //} else {
    //    frac = cw->diffuse_frac;
    //}

    // Hydraulic conductance of the entire soil-to-leaf pathway
    ktot = 1.0 / (f->total_soil_resist + 1.0 / kp);

    // Maximum transpiration rate (mmol m-2 s-1)
    emax_leaf = ktot * (s->weighted_swp - p->min_lwp);

    // Leaf transpiration (mmol m-2 s-1), i.e. ignoring boundary layer effects!
    etest = MOL_2_MMOL * (m->vpd / m->press) * cw->gsc_leaf[idx] * GSVGSC;

    // leaf water potential (MPa)
    cw->lwp_leaf[idx] = calc_lwp(f, s, ktot, etest);

    if (etest > emax_leaf) {
        stressed = TRUE;

        // Calculate gs (mol m-2 s-1) given emax_leaf
        gsv = MMOL_2_MOL * emax_leaf / (m->vpd / m->press);
        cw->gsc_leaf[idx] = gsv / GSVGSC;


        // Need to make sure transpiration solution is consistent, force
        // Tleaf to Tair as we aren't solving this
        cw->trans_leaf[idx] = emax_leaf * MMOL_2_MOL;

        cw->tleaf[idx] = m->tair;
        cw->tleaf_new = m->tair;
        cw->rnet_leaf[idx] = calc_leaf_net_rad(p, s, m->tair, m->vpd,
                                               cw->apar_leaf[idx] * PAR_2_SW);

        // Minimum leaf water potential reached so recalculate LWP (MPa)
        cw->lwp_leaf[idx] = calc_lwp(f, s, ktot, emax_leaf);

        // Re-solve An for the new gs
        photosynthesis_C3_emax(c, cw, m, p, s);

        // Need to calculate an effective beta to use in soil decomposition
        //cw->fwsoil_leaf[idx] = emax_leaf / etest;
        cw->fwsoil_leaf[idx] = exp(p->g1 * s->predawn_swp);
    } else {
        cw->fwsoil_leaf[idx] = 1.0;

    }



    return (stressed);
}

double calc_lwp(fluxes *f, state *s, double kl, double transpiration) {

    double lwp;

    lwp = s->weighted_swp - (transpiration / kl);
    if (lwp < -20.0) {
        lwp = -20.0;
    }

    return (lwp);
}

void unpack_solar_geometry(canopy_wk *cw, control *c) {

    // This geometry calculations are suprisingly intensive which is a waste
    // during spinup, so we are now doing this once and then we are just
    // accessing the 30-min value from the array position

    //calculate_solar_geometry(cw, p, doy, hod);
    //get_diffuse_frac(cw, doy, m->sw_rad);
    cw->cos_zenith = cw->cz_store[c->hour_idx];
    cw->elevation = cw->ele_store[c->hour_idx];
    cw->diffuse_frac = cw->df_store[c->hour_idx];

    return;
}
