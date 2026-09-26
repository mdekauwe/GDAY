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
    double doy, year, dummy2=0.0, relk;

    // Hydraulic conductance of the entire soil-to-leaf pathway
    // - this is only used in hydraulics, so set it to zero.
    // (mmol m–2 s–1 MPa–1)
    double ktot = 0.0;

    /* loop through the day */
    zero_carbon_day_fluxes(f);
    zero_water_day_fluxes(f);
    sunlight_hrs = 0;
    doy = ma->doy[c->hour_idx];
    year = ma->year[c->hour_idx];

    // reset plant water store to yesterday's value
    if (c->water_store) {
        // Assign plant hydraulic conductance (mmol m–2 s–1 MPa–1) from PLC
        // curve and stem water potential
        relk = calc_relative_weibull(cw->xylem_psi, p->p50, p->plc_shape);
        cw->plant_k = relk * p->kp;
    } else {
        // no cavitation when stem water storage not simulated
        cw->plant_k = p->kp;
    }

    for (hod = 0; hod < c->num_hlf_hrs; hod++) {
        unpack_met_data(c, f, ma, m, hod, dummy2);

        //if (year >= 2004.0 && year <=2005.0) {
        //    m->rain = 0.0;
        //}

        /* calculates diffuse frac from half-hourly incident radiation */
        unpack_solar_geometry(cw, c);

        /* Is the sun up? */
        if (cw->elevation > 0.0 && m->par > 20.0) {
            calculate_absorbed_radiation(cw, p, s, m->sw_rad, m->tair);
            calculate_top_of_canopy_leafn(cw, p, s);
            calc_leaf_to_canopy_scalar(cw, p, s);

            /* sunlit / shaded loop */
            for (cw->ileaf = 0; cw->ileaf < NUM_LEAVES; cw->ileaf++) {

                /* initialise values of Tleaf, Cs, dleaf at the leaf surface */
                initialise_leaf_surface(cw, m);
                iter = 0;

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
                            calculate_emax(c, cw, f, m, p, s, &ktot);
                        }

                        /* Calculate new Cs, dleaf, Tleaf */
                        solve_leaf_energy_balance(c, cw, f, m, p, s, ktot);

                    } else {
                        /*
                        ** No carbon gain, so gs = g0 ~ 0. Don't carry over
                        ** the water fluxes from a previous iteration or
                        ** timestep.
                        */
                        zero_leaf_water_fluxes(c, cw, s);
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
            for (cw->ileaf = 0; cw->ileaf < NUM_LEAVES; cw->ileaf++) {
                zero_leaf_water_fluxes(c, cw, s);
            }

            /* set tleaf to tair during the night */
            cw->tleaf[SUNLIT] = m->tair;
            cw->tleaf[SHADED] = m->tair;

            /*
            ** pre-dawn soil water potential (MPa), clearly one should link this
            ** the actual sun-rise :). Here 10 = 5 am, 10 is num_half_hr
            **/
            if (c->water_balance == HYDRAULICS && hod == 10) {
                s->predawn_swp = s->weighted_swp;

            }

        }


        scale_leaf_to_canopy(c, cw, s);
        if (c->water_balance == HYDRAULICS && hod == 24) {
            s->midday_lwp = cw->lwp_canopy;
            s->midday_xwp = cw->xylem_psi;
        }
        sum_hourly_carbon_fluxes(cw, f, p);

        // We need to remove the et_deficit which will come from the
        // plant storage from the water we need to extract from the soil.
        // We will add this back later to the transpiration output.
        if (c->water_balance == HYDRAULICS && c->water_store) {
            cw->trans_canopy -= cw->trans_deficit_canopy ;
            if (cw->trans_canopy < 0.0) {
                cw->trans_canopy = 0.0;
            }
        }

        calculate_water_balance_sub_daily(c, cw, f, m, nr, p, s, dummy,
                                          cw->trans_canopy, cw->omega_canopy,
                                          cw->rnet_canopy,
                                          cw->trans_deficit_canopy, year, doy);

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

    return;
}

void solve_leaf_energy_balance(control *c, canopy_wk *cw, fluxes *f, met *m,
                               params *p, state *s, double ktot) {
    /*
        Wrapper to solve conductances, transpiration and calculate a new
        leaf temperautre, vpd and Cs at the leaf surface.

        - The logic broadly follows MAESTRA code, with some restructuring.

        References
        ----------
        * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.

    */
    int    idx;
    double omega, transpiration, LE, Tdiff, gv, gbc, gh, trans_mmol, lai_leaf;

    idx = cw->ileaf;

    // floor avoids a zero conductance when one fraction has ~no leaf area,
    // its rnet and An are then ~0 too
    lai_leaf = MAX(cw->lai_leaf[idx], 0.001);

    penman_leaf_wrapper(m, p, s, cw->tleaf[idx], cw->rnet_leaf[idx],
                        cw->gsc_leaf[idx], lai_leaf, &transpiration, &LE, &gbc,
                        &gh, &gv, &omega);

    /* store in structure */
    cw->trans_leaf[idx] = transpiration;
    cw->omega_leaf[idx] = omega;

    /*
     * calculate new Cs, dleaf & tleaf
     */
    // Rn_iso - LE = cp Ma gh (Tleaf - Tair), gh including the radiative
    // conductance. MAESPA has Tdiff / 4 here, but that makes the converged
    // leaf-air difference 4x too small.
    Tdiff = (cw->rnet_leaf[idx] - LE) / (CP * MASS_AIR * gh);
    cw->tleaf_new = m->tair + Tdiff;
    cw->Cs = m->Ca - cw->an_leaf[idx] / gbc;
    cw->dleaf = cw->trans_leaf[idx] * m->press / gv;

    if (c->water_balance == HYDRAULICS) {
        // leaf water potential (MPa)
        trans_mmol = cw->trans_leaf[idx] * MOL_2_MMOL;
        cw->lwp_leaf[idx] = calc_lwp(f, s, ktot, trans_mmol);
    }

    return;
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

void zero_leaf_water_fluxes(control *c, canopy_wk *cw, state *s) {
    /*
        Reset the current leaf's water fluxes for when it isn't transpiring,
        i.e. at night or when An is ~0. With no flow the leaf water potential
        equilibrates with the soil.
    */
    int idx = cw->ileaf;

    cw->trans_leaf[idx] = 0.0;
    cw->omega_leaf[idx] = 0.0;

    if (c->water_balance == HYDRAULICS) {
        cw->trans_deficit_leaf[idx] = 0.0;
        cw->lwp_leaf[idx] = s->weighted_swp;
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
        // mmol m-2 s-1 to mol m-2 s-1, for consistency with transpiration
        cw->trans_deficit_canopy = (cw->trans_deficit_leaf[SUNLIT] +
                                   cw->trans_deficit_leaf[SHADED]) * MMOL_2_MOL;
    }


    return;
}

void sum_hourly_carbon_fluxes(canopy_wk *cw, fluxes *f, params *p) {

    /*
    ** GPP is gross, i.e. An + Rd, as NPP = CUE * GPP already accounts for
    ** all autotrophic respiration (including leaf Rd).
    ** umol m-2 s-1 -> gC m-2 30 min-1
    */
    f->gpp_gCm2 += (cw->an_canopy + cw->rd_canopy) * UMOL_TO_MOL * \
                    MOL_C_TO_GRAMS_C * SEC_2_HLFHR;
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

void calc_leaf_to_canopy_scalar(canopy_wk *cw, params *p, state *s) {
    /*
        Calculate scalar to transform beam/diffuse leaf Vcmax, Jmax and Rd values
        to big leaf values.

        - Insert eqn C6 & C7 into B5

        Parameters:
        ----------
        canopy_wk : structure
            various canopy values: in this case the sunlit or shaded LAI &
            cos_zenith angle.


        References:
        ----------
        * Wang and Leuning (1998) AFm, 91, 89-111; particularly the Appendix.
    */
    double kn = p->kn;

    // Parameters to scale up from single leaf to the big leaves
    cw->scalex[SUNLIT] = (1.0 - exp(-cw->kb * s->lai) * \
                                exp(-kn * s->lai)) / (cw->kb + kn);
    cw->scalex[SHADED] = (1.0 - exp(-kn * s->lai)) / kn - cw->scalex[SUNLIT];

    return;
}



void calculate_emax(control *c, canopy_wk *cw, fluxes *f, met *m, params *p,
                    state *s, double *ktot) {

    // Assumption that during the day transpiration cannot exceed a maximum
    // value, Emax (e_supply). At this point we've reached a leaf water
    // potential minimum. Once this point is reached transpiration, gs and A
    // are reclulated
    //
    // Reference:
    // * Duursma et al. 2008, Tree Physiology 28, 265–276

    double e_supply, e_demand, gsv;
    int    idx = cw->ileaf;

    // Hydraulic conductance of the entire soil-to-leaf pathway
    // (mmol m–2 s–1 MPa–1)
    *ktot = 1.0 / (f->total_soil_resist + 1.0 / cw->plant_k);

    // Maximum transpiration rate (mmol m-2 s-1)
    // Following Darcy's law which relates leaf transpiration to hydraulic
    // conductance of the soil-to-leaf pathway and leaf & soil water potentials.
    // Transpiration is limited in the perfectly isohydric case above the
    // critical threshold for embolism given by min_lwp.
    e_supply = MAX(0.0, *ktot * (s->weighted_swp - p->min_lwp));

    // Leaf transpiration (mmol m-2 s-1), i.e. ignoring boundary layer effects!
    e_demand = MOL_2_MMOL * (m->vpd / m->press) * cw->gsc_leaf[idx] * GSVGSC;

    if (e_demand > e_supply) {

        // Calculate gs (mol m-2 s-1) given supply (Emax)
        gsv = MMOL_2_MOL * e_supply / (m->vpd / m->press);
        cw->gsc_leaf[idx] = gsv / GSVGSC;

        // gs cannot be lower than minimum (cuticular conductance)
        if (cw->gsc_leaf[idx] < p->gs_min) {
            cw->gsc_leaf[idx] = p->gs_min;
            gsv = cw->gsc_leaf[idx] * GSVGSC;
        }

        // Need to calculate an effective beta to use in soil decomposition
        cw->fwsoil_leaf[idx] = e_supply / e_demand;
        //cw->fwsoil_leaf[idx] = exp(p->g1 * s->predawn_swp);

        // Re-solve An for the new gs
        photosynthesis_C3_emax(c, cw, m, p, s, cw->apar_leaf[idx],
                               cw->fwsoil_leaf[idx]);

    } else {

        cw->fwsoil_leaf[idx] = 1.0;
        gsv = cw->gsc_leaf[idx] * GSVGSC;

    }

    // Transpiration minus supply by soil/plant (emax) must be drawn from
    // plant reserve (mmol m-2 s-1). As long as there is sufficient soil water
    // this will be 0 as gsv will have been recalculated from the supply. There
    // will only be a deficit when the soil is empty and cuticular conductance
    // has taken over.
    cw->trans_deficit_leaf[idx] = MAX(0.0,
                                      (m->vpd / m->press) * gsv * MOL_2_MMOL -\
                                       e_supply);
    return;
}

double calc_lwp(fluxes *f, state *s, double ktot, double transpiration) {

    double lwp;

    if (ktot > 0.0) {
        lwp = s->weighted_swp - (transpiration / ktot);
    } else {
        lwp = s->weighted_swp;
    }

    // Set lower limit to LWP
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
