/* ============================================================================
* Calculates all within canopy C & water fluxes.
*
*
* NOTES:
*   - Should resturcure the code so that MATE is called from within the canopy
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

void canopy(control *c, fluxes *f, met *m, params *p, state *s) {
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
           press, vpd, par, tair, wind, Ca, sunlit_lai, shaded_lai, sw_rad,
           trans_canopy;
    double an_leaf[2], gsc_leaf[2], apar_leaf[2], trans_leaf[2], N0[2];
    int    hod, iter = 0, itermax = 100, i, dummy;

    /* loop through the day */
    zero_carbon_day_fluxes(f);
    zero_water_day_fluxes(f);
    for (hod = 0; hod < c->num_hlf_hrs; hod++) {

        calculate_solar_geometry(p, m->doy[c->hrly_idx], hod, &cos_zenith,
                                 &elevation);

        /* calculates diffuse frac from half-hourly incident radiation */
        par = m->par[c->hrly_idx];
        sw_rad = par * PAR_2_SW; /* SW_down [W/m2] = [J m-2 s-1] */
        diffuse_frac = get_diffuse_frac(m->doy[c->hrly_idx], cos_zenith,
                                        sw_rad);

        /* Is the sun up? */
        if (elevation > 0.0 && par > 50.0) {
            calculate_absorbed_radiation(p, s, par, diffuse_frac, elevation,
                                         cos_zenith, &(apar_leaf[0]),
                                         &sunlit_lai, &shaded_lai);

            /* Not sure if this quite makes sense for shaded bit? */
            calculate_top_of_canopy_leafn(p, s, sunlit_lai, shaded_lai,
                                          &(N0[0]));

            /* sunlit / shaded loop */
            for (i = 0; i < NUM_LEAVES; i++) {

                /* initialise values of Tleaf, Cs, dleaf at the leaf surface */
                tleaf = m->tair[c->hrly_idx];
                dleaf = m->vpd[c->hrly_idx] * KPA_2_PA;
                Cs = m->co2[c->hrly_idx];

                /* Leaf temperature stability loop */
                while (TRUE) {

                    if (c->ps_pathway == C3) {
                        photosynthesis_C3(c, p, s, N0[i], tleaf, apar_leaf[i],
                                          Cs, dleaf, &an_leaf[i], &gsc_leaf[i]);
                    } else {
                        /* Nothing implemented */
                        fprintf(stderr, "C4 photosynthesis not implemented\n");
                        exit(EXIT_FAILURE);
                    }

                    if (an_leaf[i] > 0.0) {
                        /* Calculate new Cs, dleaf, Tleaf */
                        solve_leaf_energy_balance(c, f, m, p, s, tleaf,
                                                  an_leaf[i], gsc_leaf[i],
                                                  apar_leaf[i], &Cs, &dleaf,
                                                  &tleaf_new, &trans_leaf[i]);
                    } else {
                        break;
                    }

                    if (iter >= itermax) {
                        fprintf(stderr, "No convergence in canopy loop\n");
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
        sum_hourly_carbon_fluxes(f, p, an_leaf, gsc_leaf, apar_leaf);
        trans_canopy = trans_leaf[SUNLIT] + trans_leaf[SHADED];
        calculate_water_balance(c, f, m, p, s, dummy, dummy, trans_canopy);

        c->hrly_idx++;
    } /* end of hour loop */


    return;
}



void solve_leaf_energy_balance(control *c, fluxes *f, met *m, params *p,
                               state *s, double tleaf, double an_leaf,
                               double gsc_leaf, double apar_leaf, double *Cs,
                               double *dleaf, double *tleaf_new,
                               double *transpiration) {
    /*
        Wrapper to solve conductances, transpiration and calculate a new
        leaf temperautre, vpd and Cs at the leaf surface.

        - The logic broadly follows MAESTRA code, with some restructuring.

        References
        ----------
        * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.

    */
    double LE; /* latent heat (W m-2) */
    double lambda, arg1, arg2, slope, gradn, gbhu, gbhf, gbh, gh, gbv, gsv, gv;
    double gbc, gamma, epsilon, omega, Tdiff, sensible_heat, rnet, ea, ema, Tk;
    double emissivity_atm, sw_rad, net_lw_rad;

    /* unpack the met data and get the units right */
    double press = m->press[c->hrly_idx] * KPA_2_PA;
    double vpd = m->vpd[c->hrly_idx] * KPA_2_PA;
    double tair = m->tair[c->hrly_idx];
    double wind = m->wind[c->hrly_idx];
    double Ca = m->co2[c->hrly_idx];

    /*
        extinction coefficient for diffuse radiation and black leaves
        (m2 ground m2 leaf)
    */
    double kd = 0.8;
    Tk = m->tair[c->hrly_idx] + DEG_TO_KELVIN;

    /* Radiation conductance (mol m-2 s-1) */
    gradn = calc_radiation_conductance(tair);

    /* Boundary layer conductance for heat - single sided, forced
       convection (mol m-2 s-1) */
    gbhu = calc_bdn_layer_forced_conduct(tair, press, wind, p->leaf_width);

    /* Boundary layer conductance for heat - single sided, free convection */
    gbhf = calc_bdn_layer_free_conduct(tair, tleaf, press, p->leaf_width);

    /* Total boundary layer conductance for heat */
    gbh = gbhu + gbhf;

    /* Total conductance for heat - two-sided */
    gh = 2.0 * (gbh + gradn);

    /* Total conductance for water vapour */
    gbv = GBVGBH * gbh;
    gsv = GSVGSC * gsc_leaf;
    gv = (gbv * gsv) / (gbv + gsv);
    gbc = gbh / GBHGBC;

    /* Isothermal net radiation (Leuning et al. 1995, Appendix) */
    ea = calc_sat_water_vapour_press(tair) - vpd;

    /* apparent emissivity for a hemisphere radiating at air temp eqn D4 */
    emissivity_atm = 0.642 * pow((ea / Tk), (1.0 / 7.0));
    sw_rad = apar_leaf * PAR_2_SW; /* W m-2 */

    /* isothermal net LW radiaiton at top of canopy, assuming emissivity of
       the canopy is 1 */
    net_lw_rad = (1.0 - emissivity_atm) * SIGMA * pow(Tk, 4.0);
    rnet = p->leaf_abs * sw_rad - net_lw_rad * kd * exp(-kd * s->lai);

    /* Penman-Monteith equation */
    *transpiration = penman_leaf(press, rnet, vpd, tair, gh, gv, gbv, gsv, &LE);

    /* sensible heat exchanged between leaf and surroundings */
    sensible_heat = (1.0 / (1.0 + gradn / gbh)) * (rnet - LE);

    /*
    ** calculate new dleaf, tleaf and Cs
    */
    /* Temperature difference between the leaf surface and the air */
    Tdiff = (rnet - LE) / (CP * MASS_AIR * gh);
    *tleaf_new = tair + Tdiff / 4.;
    *Cs = Ca - an_leaf / gbc;                /* CO2 conc at the leaf surface */
    *dleaf = *transpiration * press / gv;         /* VPD at the leaf surface */

    return;
}

double calc_radiation_conductance(double tair) {
    /*  Returns the 'radiation conductance' at given temperature.

        Units: mol m-2 s-1

        References:
        -----------
        * Formula from Ying-Ping's version of Maestro, cf. Wang and Leuning
          1998, Table 1,
        * See also Jones (1992) p. 108.
        * And documented in Medlyn 2007, equation A3, although I think there
          is a mistake. It should be Tk**3 not Tk**4, see W & L.
    */
    double grad;
    double Tk;

    Tk = tair + DEG_TO_KELVIN;
    grad = 4.0 * SIGMA * (Tk * Tk * Tk) * LEAF_EMISSIVITY / (CP * MASS_AIR);

    return (grad);
}

double calc_bdn_layer_forced_conduct(double tair, double press, double wind,
                                     double leaf_width) {
    /*
        Boundary layer conductance for heat - single sided, forced convection
        (mol m-2 s-1)
        See Leuning et al (1995) PC&E 18:1183-1200 Eqn E1
    */
    double cmolar, Tk, gbh;

    Tk = tair + DEG_TO_KELVIN;
    cmolar = press / (RGAS * Tk);
    gbh = 0.003 * sqrt(wind / leaf_width) * cmolar;

    return (gbh);
}

double calc_bdn_layer_free_conduct(double tair, double tleaf, double press,
                                   double leaf_width) {
    /*
        Boundary layer conductance for heat - single sided, free convection
        (mol m-2 s-1)
        See Leuning et al (1995) PC&E 18:1183-1200 Eqns E3 & E4
    */
    double cmolar, Tk, gbh, grashof, leaf_width_cubed;
    double tolerance = 1E-08;

    Tk = tair + DEG_TO_KELVIN;
    cmolar = press / (RGAS * Tk);
    leaf_width_cubed = leaf_width * leaf_width * leaf_width;

    if (float_eq((tleaf - tair), 0.0)) {
        gbh = 0.0;
    } else {
        grashof = 1.6E8 * fabs(tleaf - tair) * leaf_width_cubed;
        gbh = 0.5 * DHEAT * pow(grashof, 0.25) / leaf_width * cmolar;
    }

    return (gbh);
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




void calculate_top_of_canopy_leafn(params *p, state *s, double sunlit_lai,
                                   double shaded_lai, double *N0)  {

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
    double Ntot_sun, Ntot_sha, Ntot, N0_canopy;
    double k = p->kext;

    /* leaf mass per area (g C m-2 leaf) */
    double LMA = 1.0 / p->sla * p->cfracts * KG_AS_G;

    if (s->lai > 0.0) {

        /* the total amount of nitrogen in sunlit/shaded parts of canopy */
        Ntot = s->shootnc * LMA * s->lai;
        Ntot_sun = s->shootnc * LMA * sunlit_lai;
        Ntot_sha = s->shootnc * LMA * shaded_lai;

        /*printf("%.10lf %.10lf %.10lf\n", Ntot, Ntot_sun , Ntot_sha);*/

        /* top of canopy leaf N in the shaded/sunlit part of canopy (gN m-2) */
        N0_canopy = Ntot * k / (1.0 - exp(-k * s->lai));
        *(N0+SUNLIT) = Ntot_sun * k / (1.0 - exp(-k * sunlit_lai));
        *(N0+SHADED) = Ntot_sha * k / (1.0 - exp(-k * shaded_lai));

        /**(N0+SHADED) = N0_canopy - *(N0+SUNLIT);*/

        /*printf("%.10lf %.10lf %.10lf\n", N0_canopy, *(N0+SUNLIT) , *(N0+SHADED));
        exit(1);*/


    } else {
        *(N0+SUNLIT) = 0.0;
        *(N0+SHADED) = 0.0;
    }
    return;
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
