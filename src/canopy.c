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
    int    n_beta = 0;
    double doy, year, dummy2=0.0, sum_beta = 0.0;
    double dT, t_lo, t_hi;

    /* loop through the day */
    zero_carbon_day_fluxes(f);
    zero_water_day_fluxes(f);
    sunlight_hrs = 0;
    doy = ma->doy[c->hour_idx];
    year = ma->year[c->hour_idx];

    for (hod = 0; hod < c->num_hlf_hrs; hod++) {
        unpack_met_data(c, f, ma, m, hod, dummy2);

        //if (year >= 2004.0 && year <=2005.0) {
        //    m->rain = 0.0;
        //}

        /* calculates diffuse frac from half-hourly incident radiation */
        unpack_solar_geometry(cw, c);

        /* Is the sun up? */
        if (cw->elevation > 0.0 && m->par > 20.0) {
            calculate_absorbed_radiation(cw, p, s, m->sw_rad, m->tair,
                                         m->lwdown);
            calc_leaf_bl_forced_conduct(cw, p, s, m);
            calculate_top_of_canopy_leafn(cw, p, s);
            calc_leaf_to_canopy_scalar(cw, p, s);

            /* sunlit / shaded loop */
            for (cw->ileaf = 0; cw->ileaf < NUM_LEAVES; cw->ileaf++) {

                /* initialise values of Tleaf, Cs, dleaf at the leaf surface */
                initialise_leaf_surface(cw, m);
                iter = 0;
                t_lo = -999.9;
                t_hi = 999.9;

                /* Leaf temperature loop */
                while (TRUE) {

                    if (c->ps_pathway != C3) {
                        /* Nothing implemented */
                        fprintf(stderr, "C4 photosynthesis not implemented\n");
                        exit(EXIT_FAILURE);
                    } else if (c->water_balance == HYDRAULICS) {
                        // Sperry profit maximisation with the plant
                        // hydraulics (gs_opt)
                        gs_opt_leaf(c, cw, m, p, s);
                    } else {
                        photosynthesis_C3(c, cw, m, p, s);
                    }

                    if (cw->an_leaf[cw->ileaf] > 1E-04) {

                        /* Calculate new Cs, dleaf, Tleaf */
                        solve_leaf_energy_balance(c, cw, f, m, p, s);

                    } else {
                        /*
                        ** No carbon gain, so gs = g0 ~ 0. Don't carry over
                        ** the water fluxes from a previous iteration or
                        ** timestep.
                        */
                        zero_leaf_water_fluxes(c, cw, p, s);
                        break;
                    }

                    if (iter >= itermax) {
                        fprintf(stderr, "No convergence in canopy loop: "
                                "%.0f doy %.0f hod %d leaf %d Tleaf %.2f "
                                "new %.2f Tair %.2f\n", year, doy, hod,
                                cw->ileaf, cw->tleaf[cw->ileaf],
                                cw->tleaf_new, m->tair);
                        exit(EXIT_FAILURE);
                    }
                    dT = cw->tleaf_new - cw->tleaf[cw->ileaf];
                    if (fabs(dT) < 0.02) {
                        break;
                    }

                    /*
                    ** Update temperature & do another iteration. Each
                    ** iteration tells us which side of the solution Tleaf is
                    ** on (the energy balance wants it warmer or cooler), so
                    ** keep a bracket. Until both sides are known, under-relax
                    ** (move half way to the energy balance value), which
                    ** damps the oscillation for large leaves / low wind; then
                    ** bisect. The bisection also copes with the energy
                    ** balance jumping between gs_opt Ci grid points either
                    ** side of the solution.
                    */
                    if (dT > 0.0) {
                        t_lo = MAX(t_lo, cw->tleaf[cw->ileaf]);
                    } else {
                        t_hi = MIN(t_hi, cw->tleaf[cw->ileaf]);
                    }
                    if (t_lo > -900.0 && t_hi < 900.0) {
                        if (t_hi - t_lo < 0.01) {
                            break;
                        }
                        cw->tleaf[cw->ileaf] = 0.5 * (t_lo + t_hi);
                    } else {
                        cw->tleaf[cw->ileaf] += 0.5 * dT;
                    }
                    iter++;
                } /* end of leaf temperature loop */

                /* the gs_opt water stress factor, if the leaf transpired */
                if (c->water_balance == HYDRAULICS &&
                    cw->an_leaf[cw->ileaf] > 1E-04) {
                    cw->fwsoil_leaf[cw->ileaf] = gs_opt_beta(c, cw, m, p, s);
                    sum_beta += cw->fwsoil_leaf[cw->ileaf];
                    n_beta++;
                }

            } /* end of sunlit/shaded leaf loop */

        } else {

            zero_hourly_fluxes(cw);
            calculate_soil_net_radiation_night(cw, p, s, m->tair, m->lwdown);
            for (cw->ileaf = 0; cw->ileaf < NUM_LEAVES; cw->ileaf++) {
                zero_leaf_water_fluxes(c, cw, p, s);
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
            s->midday_plc = p->kp > 0.0 ?
                            100.0 * (1.0 - cw->kl_canopy / p->kp) : 0.0;
        }
        sum_hourly_carbon_fluxes(cw, f, p);

        calculate_water_balance_sub_daily(c, cw, f, m, nr, p, s, dummy,
                                          cw->trans_canopy, cw->omega_canopy,
                                          cw->rnet_canopy, year, doy);

        if (c->print_options == SUBDAILY && c->spin_up == FALSE) {
            write_subdaily_outputs_ascii(c, cw, s, year, doy, hod);
        }
        c->hour_idx++;
        if (cw->elevation > 0.0 && m->par > 20.0) {
            sunlight_hrs++;
        }
    } /* end of hour loop */

    /* work out average omega for the day over sunlight hours */
    if (sunlight_hrs > 0) {
        f->omega /= sunlight_hrs;
    }

    if (c->water_stress && c->water_balance == HYDRAULICS) {
        // Daytime mean of the gs_opt stress factor (A / A with wet soil),
        // used in soil decomposition and allocation. Unchanged if nothing
        // transpired.
        if (n_beta > 0) {
            s->wtfac_root = sum_beta / (double)n_beta;
            s->wtfac_topsoil = s->wtfac_root;
        }
    } else if (c->water_stress) {
        // Calculate the soil moisture availability factors [0,1] in the
        // topsoil and the entire root zone
        calculate_soil_water_fac(c, p, s);
    } else {
        s->wtfac_topsoil = 1.0;
        s->wtfac_root = 1.0;
    }

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
    double omega, transpiration, LE, Tdiff, gv, gbc, gh, trans_mmol, lai_leaf;

    idx = cw->ileaf;

    lai_leaf = cw->lai_leaf[idx];

    penman_leaf_wrapper(m, p, s, cw->tleaf[idx], cw->rnet_leaf[idx],
                        cw->gsc_leaf[idx], cw->gbhu[idx], cw->gradis[idx],
                        lai_leaf, &transpiration, &LE, &gbc, &gh, &gv, &omega);

    /* store in structure */
    cw->trans_leaf[idx] = transpiration;
    cw->omega_leaf[idx] = omega;

    /*
     * calculate new Cs, dleaf & tleaf
     */
    // Rn_iso - LE = cp Ma gh (Tleaf - Tair), gh including the (big-leaf)
    // radiative conductance. MAESPA has Tdiff / 4 here, but that makes the converged
    // leaf-air difference 4x too small.
    Tdiff = (cw->rnet_leaf[idx] - LE) / (CP * MASS_AIR * gh);
    cw->tleaf_new = m->tair + Tdiff;
    cw->Cs = m->Ca - cw->an_leaf[idx] / gbc;
    cw->dleaf = cw->trans_leaf[idx] * m->press / gv;

    if (c->water_balance == HYDRAULICS) {
        // leaf water potential (MPa) supplying the transpiration
        trans_mmol = cw->trans_leaf[idx] * MOL_2_MMOL;
        cw->lwp_leaf[idx] = gs_opt_psi_leaf(cw, p, s, trans_mmol,
                                            &cw->kl_leaf[idx]);
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

void zero_leaf_water_fluxes(control *c, canopy_wk *cw, params *p, state *s) {
    /*
        Reset the current leaf's water fluxes for when it isn't transpiring,
        i.e. at night or when An is ~0. With no flow the leaf water potential
        equilibrates with the soil.
    */
    int idx = cw->ileaf;

    cw->trans_leaf[idx] = 0.0;
    cw->omega_leaf[idx] = 0.0;

    if (c->water_balance == HYDRAULICS) {
        // no flow, the leaf is at the root zone water potential
        cw->lwp_leaf[idx] = gs_opt_psi_leaf(cw, p, s, 0.0, &cw->kl_leaf[idx]);
    }

    return;
}

void scale_leaf_to_canopy(control *c, canopy_wk *cw, state *s) {

    cw->an_canopy = cw->an_leaf[SUNLIT] + cw->an_leaf[SHADED];
    cw->rd_canopy = cw->rd_leaf[SUNLIT] + cw->rd_leaf[SHADED];
    cw->gsc_canopy = cw->gsc_leaf[SUNLIT] + cw->gsc_leaf[SHADED];
    cw->apar_canopy = cw->apar_leaf[SUNLIT] + cw->apar_leaf[SHADED];
    cw->trans_canopy = cw->trans_leaf[SUNLIT] + cw->trans_leaf[SHADED];
    cw->omega_canopy = (cw->omega_leaf[SUNLIT] + cw->omega_leaf[SHADED]) / 2.0;
    cw->rnet_canopy = cw->rnet_leaf[SUNLIT] + cw->rnet_leaf[SHADED];

    if (c->water_balance == HYDRAULICS) {
        // leaf area weighted means
        double lai = cw->lai_leaf[SUNLIT] + cw->lai_leaf[SHADED];
        if (lai > 0.0) {
            cw->lwp_canopy = (cw->lwp_leaf[SUNLIT] * cw->lai_leaf[SUNLIT] +
                              cw->lwp_leaf[SHADED] * cw->lai_leaf[SHADED]) /
                             lai;
            cw->kl_canopy = (cw->kl_leaf[SUNLIT] * cw->lai_leaf[SUNLIT] +
                             cw->kl_leaf[SHADED] * cw->lai_leaf[SHADED]) / lai;
        } else {
            cw->lwp_canopy = s->weighted_swp;
            cw->kl_canopy = cw->kl_leaf[SUNLIT];
        }
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
