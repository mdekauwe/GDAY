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
#include "water_balance.h"

void canopy(control *c, fluxes *f, met *m, params *p, state *s,
            double *day_tsoil, double *day_ndep) {
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
    double Cs, dleaf, tleaf, tleaf_new, trans_hlf_hr, leafn, fc, cos_zenith,
           elevation, direct_apar, diffuse_apar, diffuse_frac, rnet=0.0,
           press, vpd, par, tair, wind, Ca, sw_rad, N0, rnet_canopy,
           trans_canopy, omega_canopy;
    double an_leaf[2], gsc_leaf[2], apar_leaf[2], trans_leaf[2], rnet_leaf[2],
           sunlit_shaded_lai[2], omega_leaf[2];
    int    hod, iter = 0, itermax = 100, i, dummy, sunlight_hrs;

    /* loop through the day */
    zero_carbon_day_fluxes(f);
    zero_water_day_fluxes(f);
    sunlight_hrs = 0;
    *day_tsoil = 0.0;
    *day_ndep = 0.0;
    for (hod = 0; hod < c->num_hlf_hrs; hod++) {

        *day_tsoil += m->tsoil[c->hrly_idx];
        *day_ndep += m->ndep[c->hrly_idx];

        calculate_solar_geometry(p, m->doy[c->hrly_idx], hod, &cos_zenith,
                                 &elevation);

        /* calculates diffuse frac from half-hourly incident radiation */
        par = m->par[c->hrly_idx];
        sw_rad = par * PAR_2_SW; /* SW_down [W/m2] = [J m-2 s-1] */
        diffuse_frac = get_diffuse_frac(m->doy[c->hrly_idx], cos_zenith,
                                        sw_rad);

        /* Is the sun up? */
        if (elevation > 0.0 && par > 20.0) {
            calculate_absorbed_radiation(p, s, par, diffuse_frac, elevation,
                                         cos_zenith, &(apar_leaf[0]),
                                         &(sunlit_shaded_lai[0]));

            /* Not sure if this quite makes sense for shaded bit? */
            N0 = calculate_top_of_canopy_leafn(p, s);

            /* sunlit / shaded loop */
            for (i = 0; i < NUM_LEAVES; i++) {

                /* initialise values of Tleaf, Cs, dleaf at the leaf surface */
                tleaf = m->tair[c->hrly_idx];
                dleaf = m->vpd[c->hrly_idx] * KPA_2_PA;
                Cs = m->co2[c->hrly_idx];

                /* Leaf temperature stability loop */
                while (TRUE) {

                    if (c->ps_pathway == C3) {
                        photosynthesis_C3(c, p, s, N0, tleaf, apar_leaf[i],
                                          Cs, dleaf, &an_leaf[i], &gsc_leaf[i],
                                          sunlit_shaded_lai[i], i, cos_zenith);
                    } else {
                        /* Nothing implemented */
                        fprintf(stderr, "C4 photosynthesis not implemented\n");
                        exit(EXIT_FAILURE);
                    }

                    if (an_leaf[i] > 1E-04) {
                        /* Calculate new Cs, dleaf, Tleaf */
                        solve_leaf_energy_balance(c, f, m, p, s, tleaf,
                                                  an_leaf[i], gsc_leaf[i],
                                                  apar_leaf[i], &Cs, &dleaf,
                                                  &tleaf_new, &trans_leaf[i],
                                                  &omega_leaf[i],
                                                  &rnet_leaf[i]);
                    } else {
                        break;
                    }

                    if (iter >= itermax) {
                        fprintf(stderr, "No convergence in canopy loop:\n");
                        exit(EXIT_FAILURE);
                    } else if (fabs(tleaf - tleaf_new) < 0.02) {
                        break;
                    }

                    /* Update temperature & do another iteration */
                    tleaf = tleaf_new;
                    iter++;

                } /* end of leaf temperature stability loop */

            } /* end of sunlit/shaded leaf loop */
        } else {
            zero_hourly_fluxes(&(an_leaf[0]), &(gsc_leaf[0]), &(trans_leaf[0]));
        }
        /*
        an_leaf[SUNLIT] *= sunlit_shaded_lai[SUNLIT];
        an_leaf[SHADED] *= sunlit_shaded_lai[SHADED];
        trans_canopy = (trans_leaf[SUNLIT]*sunlit_shaded_lai[SUNLIT]) + (trans_leaf[SHADED]*sunlit_shaded_lai[SHADED]);
        */
        sum_hourly_carbon_fluxes(f, p, an_leaf, gsc_leaf, apar_leaf);
        trans_canopy = trans_leaf[SUNLIT] + trans_leaf[SHADED];
        omega_canopy = (omega_leaf[SUNLIT] + omega_leaf[SHADED]) / 2.0;
        rnet_canopy = rnet_leaf[SUNLIT] + rnet_leaf[SHADED];
        calculate_water_balance(c, f, m, p, s, dummy, dummy, trans_canopy,
                                omega_canopy, rnet_canopy);

        c->hrly_idx++;
        sunlight_hrs++;
    } /* end of hour loop */

    /* work out average omega for the day over sunlight hours */
    f->omega /= sunlight_hrs;

    /* work out average omega for the day, including the night */
    *day_tsoil /= c->num_hlf_hrs;

    if (c->water_stress) {
        /* Calculate the soil moisture availability factors [0,1] in the
           topsoil and the entire root zone */
        calculate_soil_water_fac(c, p, s);
    } else {
        /* really this should only be a debugging option! */
        s->wtfac_topsoil = 1.0;
        s->wtfac_root = 1.0;
    }

    return;
}



void solve_leaf_energy_balance(control *c, fluxes *f, met *m, params *p,
                               state *s, double tleaf, double an_leaf,
                               double gsc_leaf, double apar_leaf, double *Cs,
                               double *dleaf, double *tleaf_new,
                               double *transpiration, double *omega,
                               double *rnet) {
    /*
        Wrapper to solve conductances, transpiration and calculate a new
        leaf temperautre, vpd and Cs at the leaf surface.

        - The logic broadly follows MAESTRA code, with some restructuring.

        References
        ----------
        * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.

    */
    double LE; /* latent heat (W m-2) */
    double Tdiff, press, vpd, tair, wind, Ca, gv, gbc, gh, Tk, sw_rad;

    /* unpack the met data and get the units right */
    press = m->press[c->hrly_idx] * KPA_2_PA;
    vpd = m->vpd[c->hrly_idx] * KPA_2_PA;
    tair = m->tair[c->hrly_idx];
    wind = m->wind[c->hrly_idx];
    Ca = m->co2[c->hrly_idx];
    sw_rad = apar_leaf * PAR_2_SW; /* W m-2 */

    *rnet = calc_leaf_net_rad(p, s, tair, vpd, sw_rad);
    penman_leaf_wrapper(p, s, press, vpd, tair, tleaf, wind, *rnet, gsc_leaf,
                        transpiration, &LE, &gbc, &gh, &gv, omega);

    /*
    ** calculate new dleaf, tleaf and Cs
    */
    /* Temperature difference between the leaf surface and the air */
    Tdiff = (*rnet - LE) / (CP * MASS_AIR * gh);
    *tleaf_new = tair + Tdiff / 4.;
    *Cs = Ca - an_leaf / gbc;                /* CO2 conc at the leaf surface */
    *dleaf = *transpiration * press / gv;         /* VPD at the leaf surface */

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




double calculate_top_of_canopy_leafn(params *p, state *s) {

    /*
    Calculate the N at the top of the canopy (g N m-2), N0.

    Returns:
    -------
    N0 : float (gN m-2)
        leaf nitrogen content at the top of the canopy

    References:
    -----------
    * Chen et al 93, Oecologia, 93,63-69.

    */
    double Ntot, N0;
    double kn = 0.3; /* extinction coefficent for Nitrogen less steep */

    /* leaf mass per area (g C m-2 leaf) */
    double LMA = 1.0 / p->sla * p->cfracts * KG_AS_G;

    if (s->lai > 0.0) {
        /* the total amount of nitrogen in the canopy */
        Ntot = s->shootnc * LMA * s->lai;

        /* top of canopy leaf N (gN m-2) */
        N0 = Ntot * kn / (1.0 - exp(-kn * s->lai));
        /*printf("%lf %lf %lf %f\n", s->lai, s->shootnc, Ntot, N0);*/

    } else {
        N0 = 0.0;
    }

    return (N0);
}

void zero_hourly_fluxes(double *an_leaf, double *gsc_leaf,
                        double *trans_leaf) {

    int i;

    /* sunlit / shaded loop */
    for (i = 0; i < NUM_LEAVES; i++) {
        *(an_leaf+i) = 0.0;
        *(gsc_leaf+i) = 0.0;
        *(trans_leaf+i) = 0.0;
    }

    return;
}

void sum_hourly_carbon_fluxes(fluxes *f, params *p, double *an_leaf,
                              double *gsc_leaf, double *apar_leaf) {

    double an_canopy = 0.0, gsc_canopy = 0.0, apar_canopy = 0.0;;
    int    i;

    /* scale to canopy */
    for (i = 0; i < NUM_LEAVES; i++) {
        an_canopy += *(an_leaf+i);
        gsc_canopy += *(gsc_leaf+i);
        apar_canopy += *(apar_leaf+i);
    }

    /* umol m-2 s-1 -> gC m-2 30 min-1 */
    f->gpp_gCm2 += an_canopy * UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HLFHR;
    f->npp_gCm2 = f->gpp_gCm2 * p->cue;
    f->gpp = f->gpp_gCm2 * GRAM_C_2_TONNES_HA;
    f->npp = f->npp_gCm2 * GRAM_C_2_TONNES_HA;
    f->auto_resp = f->gpp - f->npp;
    f->apar += apar_canopy;
    f->gs_mol_m2_sec += gsc_canopy;

    return;
}
