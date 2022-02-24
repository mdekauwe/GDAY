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
    int    debug = TRUE, k;
    double doy, year, dummy2=0.0, previous_sw, current_sw, gsv;
    double previous_cs, current_cs, relk;
    double Anet_opt[c->resolution], gsc_opt[c->resolution];
    // Hydraulic conductance of the entire soil-to-leaf pathway
    // - this is only used in hydraulics, so set it to zero.
    // (mmol m–2 s–1 MPa–1)
    double ktot = 0.0;

    /* loop through the day */
    zero_carbon_day_fluxes(f);
    zero_water_day_fluxes(f);
    previous_sw = s->pawater_topsoil + s->pawater_root;
    previous_cs = s->canopy_store;
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

                /* Leaf temperature loop */
                while (TRUE) {

                    if (c->ps_pathway == C3 && c->water_balance != GS_OPT) {
                        photosynthesis_C3(c, cw, m, p, s);
                    } else if (c->ps_pathway == C3 && c->water_balance == GS_OPT) {
                        photosynthesis_C3_opt(c, cw, m, p, s, Anet_opt, gsc_opt);

                        //for (k=0; k<c->resolution; k++) {
                        //    printf("%f %f\n", Anet_opt[k], gsc_opt[k]);
                        //}

                    } else {
                        /* Nothing implemented */
                        fprintf(stderr, "C4 photosynthesis not implemented\n");
                        exit(EXIT_FAILURE);
                    }

                    if (c->water_balance == GS_OPT) {
                        calculate_gs_E_psi_leaf(c, cw, f, m, p, s, &ktot,
                                                Anet_opt, gsc_opt);

                        if (cw->an_leaf[cw->ileaf] > 1E-04) {
                            /* Calculate new Cs, dleaf, Tleaf */
                            solve_leaf_energy_balance(c, cw, f, m, p, s, ktot);
                            //printf("**%f %f\n\n", cw->an_leaf[cw->ileaf], cw->gsc_leaf[cw->ileaf]);
                        } else {
                            break;
                        }


                    } else {


                        if (cw->an_leaf[cw->ileaf] > 1E-04) {

                            if (c->water_balance == HYDRAULICS) {
                                // Ensure transpiration does not exceed Emax, if it
                                // does we recalculate gs and An
                                calculate_emax(c, cw, f, m, p, s, &ktot);
                            }

                            /* Calculate new Cs, dleaf, Tleaf */
                            solve_leaf_energy_balance(c, cw, f, m, p, s, ktot);

                        } else {
                            break;
                        }
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

    //current_sw = s->pawater_topsoil + s->pawater_root;
    //current_cs = s->canopy_store;
    //check_water_balance(c, f, s, previous_sw, current_sw, previous_cs,
    //                    current_cs, year, doy);

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
    double omega, transpiration, LE, Tdiff, gv, gbc, gh, sw_rad, trans_mmol;

    if (c->water_balance == HYDRAULICS) {
        idx = cw->ileaf;
        penman_leaf_wrapper(m, p, s, cw->tleaf[idx], cw->rnet_leaf[idx],
                            cw->gsc_leaf[idx], &transpiration, &LE, &gbc, &gh,
                            &gv, &omega);

        /* store in structure */
        cw->trans_leaf[idx] = transpiration;
        cw->omega_leaf[idx] = omega;
    }

    /*
     * calculate new Cs, dleaf & tleaf
     */
    Tdiff = (cw->rnet_leaf[idx] - LE) / (CP * MASS_AIR * gh);
    cw->tleaf_new = m->tair + Tdiff / 4.0;
    cw->Cs = m->Ca - cw->an_leaf[idx] / gbc;
    cw->dleaf = cw->trans_leaf[idx] * m->press / gv;

    if (c->water_balance == HYDRAULICS) {
        // leaf water potential (MPa)
        trans_mmol = cw->trans_leaf[idx] * MOL_2_MMOL;
        cw->lwp_leaf[idx] = calc_lwp(f, s, ktot, trans_mmol);
    }

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
        // mmol m-2 s-1 to mol m-2 s-1, for consistency with transpiration
        cw->trans_deficit_canopy = (cw->trans_deficit_leaf[SUNLIT] +
                                   cw->trans_deficit_leaf[SHADED]) * MMOL_2_MOL;
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

void calculate_gs_E_psi_leaf(control *c, canopy_wk *cw, fluxes *f, met *m,
                             params *p, state *s, double *ktot, double *An,
                             double *gsc) {


    double e_supply, e_demand, gsv, frac;
    int    leaf_idx = cw->ileaf;
    int    k, N, idx=0;
    double e_leaf[c->resolution], psi_leaf[c->resolution], Kc[c->resolution];
    double gain[c->resolution], cost[c->resolution], profit[c->resolution];
    double RGSWC = 1.57;
    double weibull, max_an, max_profit, p_crit, kcmax;
    N = c->resolution;

    // Hydraulic conductance of the entire soil-to-leaf pathway
    // (mmol m–2 s–1 MPa–1)
    //*ktot = 1.0 / (f->total_soil_resist + 1.0 / cw->plant_k);

    // Canopy xylem pressure (P_crit) MPa, beyond which tree
    // desiccates (Ecrit), MPa
    p_crit = -p->b_plant * pow(log(cw->plant_k / p->Kcrit), (1.0 / p->c_plant));

    // Plant hydraulic conductance (mmol m-2 leaf s-1 MPa-1)
    weibull = exp(-pow((-s->weighted_swp / p->b_plant), p->c_plant));
    weibull = MIN(1.0, MAX(1.0E-09, weibull));
    kcmax = cw->plant_k * weibull;

    // Find Max An, we need this to normalise by below
    max_an = -9999.;
    for (k=0; k<N; k++) {

        // Ensure we don't check for profit in bad psi_leaf search space
        if (psi_leaf[k] <= s->weighted_swp || psi_leaf[k] >= p_crit) {

            if (An[k] > max_an) {
                max_an = An[k];
            }
        }
    }

    if (max_an > -9000.) {


        max_profit = -9999.;
        for (k=0; k<N; k++) {

            // Assuming perfect coupling, infer E_sun/sha from gsc. NB. as we're
            // iterating, Tleaf will change and so VPD, maintaining energy
            // balance
            e_leaf[k] = gsc[k] * RGSWC / m->press * cw->dleaf; // mol H2O m-2 s-1


            // Rescale from canopy to leaf..as e_leaf is E_sun/sha i.e. big-leaf to
            // unit leaf, mmol m-2 s-1
            e_leaf[k] *= MOL_2_MMOL / cw->scalex[leaf_idx];

            // Infer the matching leaf water potential (MPa).
            psi_leaf[k] = s->weighted_swp - e_leaf[k] / cw->plant_k;
            //psi_leaf[k] = MAX(-20.0, psi_leaf[k]);

            // Soil–plant hydraulic conductance at canopy xylem pressure
            // mmol m-2 s-1 MPa-1
            weibull = exp(-pow((-psi_leaf[k] / p->b_plant), p->c_plant));
            weibull = MIN(1.0, MAX(1.0E-09, weibull));
            Kc[k] = cw->plant_k * weibull;

            // normalised gain (-)
            gain[k] = An[k] / max_an;
            gain[k] = MIN(1.0, MAX(1.0E-09, gain[k]));

            // normalised cost (-)
            cost[k] = (kcmax - Kc[k]) / (kcmax - p->Kcrit);
            cost[k] = MIN(1.0, MAX(1.0E-09, cost[k]));

            // Locate maximum profit
            profit[k] = gain[k] - cost[k];


            //printf("* %f %f %f\n", profit[k], gain[k], cost[k]);


            // Ensure we don't check for profit in bad psi_leaf search space
            if (psi_leaf[k] <= s->weighted_swp || psi_leaf[k] >= p_crit) {

                if (profit[k] > max_profit) {
                    max_profit = profit[k];
                    idx = k;

                }
            }

            //printf("** %f %f %d %f %f %f\n", max_profit, profit[k], idx, psi_leaf[k], s->weighted_swp, p_crit);


        }



        // Save optimised values
        cw->an_leaf[leaf_idx] = An[idx]; // umol m-2 s-1
        cw->trans_leaf[leaf_idx] = e_leaf[idx] * MMOL_2_MOL; // ! mol H2O m-2 s-1
        cw->lwp_leaf[leaf_idx] = psi_leaf[idx]; // MPa

        //printf("%d %f %f %f %f %f %f \n", idx, An[idx], e_leaf[idx], gsc[idx], psi_leaf[idx],  gain[idx], cost[idx]);
        //if (e_leaf[idx] > 87942) {
        //    for (k=0; k<N; k++) {
        //        printf("%d: %f %f %f %f\n", k, profit[k], gain[k], cost[k], max_profit);
        //    }
        //    printf("%d\n", idx);
        //    exit(1);
        //}
        //if (e_leaf[idx] > 87942) {
        //    exit(1);
        //}

        //printf("%d %f %f %f\n", idx, cw->an_leaf[leaf_idx], cw->trans_leaf[leaf_idx], cw->lwp_leaf[leaf_idx]);


    } else {
        cw->an_leaf[leaf_idx] = 0.0; // umol m-2 s-1
        cw->trans_leaf[leaf_idx] = 0.0; // ! mol H2O m-2 s-1
        cw->lwp_leaf[leaf_idx] = s->weighted_swp; // MPa
        //printf("%f %d\n", cw->an_leaf[leaf_idx], idx);
    }

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

    double e_supply, e_demand, gsv, frac;
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
