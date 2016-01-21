/* ============================================================================
* Calculates all within canopy C & water fluxes.
*
*
* NOTES:
*
*
* AUTHOR:
*   Martin De Kauwe
*
* DATE:
*   17.02.2015
*
* =========================================================================== */
#include "canopy.h"
#include "water_balance.h"

void canopy(control *c, fluxes *f, met *m, params *p, state *s,
            int project_day) {
    /*
        Two-leaf canopy module consists of two parts:
        (1) a radiation submodel to calculate PAR, NIR and thermal radiation
        (2) a coupled model of stomatal conductance, photosynthesis and
            partitioning between absorbed net radiation into sensible and
            latent heat.

        The coupled model has two leaves: sunlit & shaded under the assumption
        that the sunlit and shaded leaf is representative of all respective
        sunlit or shaded leaves within a canopy. Clearly for dense canopies this
        assumption will not hold due, but as fluxes from canopy elements at the
        base of the canopy are small, it is likely to be an acceptable error.

        References
        ----------
        * Wang & Leuning (1998) Agricultural & Forest Meterorology, 91, 89-111.
    */

    double Cs, dleaf, tleaf, new_tleaf, trans_hlf_hr, leafn, fc, ncontent,
           cos_zenith, elevation, anleaf[2], gsc, apar[2],
           direct_apar, diffuse_apar, total_rnet, diffuse_frac, rnet,
           press, vpd, par, tair, wind, Ca, sunlit_frac, acanopy;
    int    hod, iter = 0, itermax = 100, ileaf;
    long   offset;


    if (s->lai > 0.0) {
        /* average leaf nitrogen content (g N m-2 leaf) */
        leafn = (s->shootnc * p->cfracts / p->sla * KG_AS_G);

        /* total nitrogen content of the canopy */
        ncontent = leafn * s->lai;

    } else {
        ncontent = 0.0;
    }

    /* When canopy is not closed, canopy light interception is reduced
        - calculate the fractional ground cover */
    if (s->lai < p->lai_closed) {
        /* discontinuous canopies */
        fc = s->lai / p->lai_closed;
    } else {
        fc = 1.0;
    }

    /* fIPAR - the fraction of intercepted PAR = IPAR/PAR incident at the
       top of the canopy, accounting for partial closure based on Jackson
       and Palmer (1979). */
    if (s->lai > 0.0)
        s->fipar = ((1.0 - exp(-p->kext * s->lai / fc)) * fc);
    else
        s->fipar = 0.0;

    if (c->water_stress) {
        /* Calculate the soil moisture availability factors [0,1] in the
           topsoil and the entire root zone */
        calculate_soil_water_fac(c, p, s);
    } else {
        /* really this should only be a debugging option! */
        s->wtfac_topsoil = 1.0;
        s->wtfac_root = 1.0;
    }

    zero_carbon_day_fluxes(f);
    zero_water_day_fluxes(f);

    for (hod = 0; hod < c->num_hlf_hrs; hod++) {
        offset = project_day * c->num_days + hod;

        /*
            initialise values of leaf temp, leaf surface CO2 and VPD from
            air space
        */
        tleaf = m->tair[offset];                        /* Leaf temperature */
        dleaf = m->vpd[offset] * KPA_2_PA;        /* VPD at the leaf suface */
        Cs = m->co2[offset];                /* CO2 conc. at the leaf suface */
        par = m->par[offset];

        calculate_zenith_angle(p, m->doy[project_day], hod, &cos_zenith,
                               &elevation);

        /* calculates diffuse frac from half-hourly incident radiation */
        diffuse_frac = get_diffuse_frac(m->doy[offset], cos_zenith, par);

        /* Is the sun up? */
        if (elevation > 0.0) {

            calculate_absorbed_radiation(p, s, par, diffuse_frac,
                                         cos_zenith, &(apar[0]));

            /* sunlit, shaded loop */
            acanopy = 0.0;
            total_rnet = 0.0;
            for (ileaf = 0; ileaf <= 1; ileaf++) {

                while (TRUE) {

                    if (c->ps_pathway == C3) {
                        /* ********************************************
                            NEED TO PASS THE LEAF PAR NOT TOTAL PAR !!!
                            ********************************************
                        */
                        photosynthesis_C3(c, p, s, ncontent, tleaf, apar[ileaf],
                                          Cs, dleaf, &gsc, &anleaf[ileaf]);
                    } else {
                        /* Nothing implemented */
                        fprintf(stderr, "C4 photosynthesis not implemented\n");
                        exit(EXIT_FAILURE);
                    }

                    if (anleaf[ileaf] > 0.0) {
                        /* Calculate new Cs, dleaf, Tleaf */
                        solve_leaf_energy_balance(f, m, p, s, offset, tleaf,
                                                  gsc, anleaf[ileaf], apar[ileaf], &Cs, &dleaf,
                                                  &new_tleaf, &trans_hlf_hr);
                    } else {
                        trans_hlf_hr = 0.0;
                    }

                    if (iter >= itermax) {
                        fprintf(stderr, "No convergence in canopy loop\n");
                        exit(EXIT_FAILURE);
                    }

                    /* stopping criteria */
                    if (fabs(tleaf - new_tleaf) < 0.02) {
                        break;
                    }

                    /* Update temperature & do another iteration */
                    tleaf = new_tleaf;
                    iter++;

                }
                total_rnet += rnet;


            }
            /* scale leaf flux to canopy */
            sunlit_frac = (1.0 - exp(-p->kext * s->lai)) / p->kext;
            acanopy = sunlit_frac * anleaf[SUNLIT];
            acanopy += (s->lai - sunlit_frac) * anleaf[SHADED];

            update_daily_carbon_fluxes(f, p, acanopy);
            /* transpiration mol m-2 s-1 to mm 30 min-1 */
            trans_hlf_hr *= MOLE_WATER_2_G_WATER * G_TO_KG * SEC_2_HLFHR;
            calculate_sub_daily_water_balance(c, f, m, p, s, offset,
                                              trans_hlf_hr, total_rnet);

            /*printf("* %lf %lf %lf %lf\n", m->par[offset], acanopy, anleaf[SUNLIT], anleaf[SHADED]);*/





        } else {
            /* set time slot photosynthesis/respiration to be zero, but we
               still need to calc the full water balance */

        }

    }
    exit(1);

    return;

}

void calculate_absorbed_radiation(params *p, state *s, double par,
                                  double diffuse_frac, double cos_zenith,
                                  double *apar) {
    /*
        Calculate absorded direct (beam) and diffuse radiation
    */

    /*  Calculate diffuse radiation absorbed directly. */
    apar[SHADED] = par * diffuse_frac * (1.0 - exp(-p->kext * s->lai));

    /* Calculate beam radiation absorbed by sunlit leaf area. */
    apar[SUNLIT] = par * (1.0 - diffuse_frac) / cos_zenith * p->leaf_abs;
    apar[SUNLIT] += apar[SHADED];

    printf("%lf %lf %lf %lf\n", par, apar[SUNLIT], apar[SHADED], apar[SUNLIT]+ apar[SHADED]);
    return;
}

void solve_leaf_energy_balance(fluxes *f, met *m, params *p, state *s,
                               long offset, double tleaf, double gsc,
                               double anleaf, double apar, double *Cs,
                               double *dleaf, double *new_tleaf, double *et) {
    /*
        Coupled model wrapper to solve photosynthesis, stomtal conductance
        and radiation paritioning.

    */
    double LE; /* latent heat (W m-2) */
    double Rspecifc_dry_air = 287.058; /* J kg-1 K-1 */
    double lambda, arg1, arg2, slope, gradn, gbhu, gbhf, gbh, gh, gbv, gsv, gv;
    double gbc, gamma, epsilon, omega, Tdiff, sensible_heat;

    /* unpack the met data and get the units right */
    double press = m->press[offset] * KPA_2_PA;
    double vpd = m->vpd[offset] * KPA_2_PA;
    double tair = m->tair[offset];
    double wind = m->wind[offset];
    double Ca = m->co2[offset];


    double rnet, ea, ema, Tk, sw_rad;
    double emissivity_atm;


    /*
        extinction coefficient for diffuse radiation and black leaves
        (m2 ground m2 leaf)
    */
    double kd = 0.8, net_lw_rad;

    Tk = m->tair[offset] + DEG_TO_KELVIN;

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
    gsv = GSVGSC * gsc;
    gv = (gbv * gsv) / (gbv + gsv);

    gbc = gbh / GBHGBC;

    /* Isothermal net radiation (Leuning et al. 1995, Appendix) */
    ea = calc_sat_water_vapour_press(tair) - vpd;

    /* apparent emissivity for a hemisphere radiating at air temp eqn D4 */
    emissivity_atm = 0.642 * pow((ea / Tk), (1.0 / 7.0));
    sw_rad = apar * PAR_2_SW; /* W m-2 */

    /* isothermal net LW radiaiton at top of canopy, assuming emissivity of
       the canopy is 1 */
    net_lw_rad = (1.0 - emissivity_atm) * SIGMA * pow(Tk, 4.0);
    rnet = p->leaf_abs * sw_rad - net_lw_rad * kd * exp(-kd * s->lai);

    /* Penman-Monteith equation */
    *et = penman_leaf(press, rnet, vpd, tair, gh, gv, gbv, gsv, &LE);

    /* sensible heat exchanged between leaf and surroundings */
    sensible_heat = (1.0 / (1.0 + gradn / gbh)) * (rnet - LE);

    /*
    ** calculate new dleaf, tleaf and Cs
    */

    /* Temperature difference between the leaf surface and the air */
    Tdiff = (rnet - LE) / (CP * MASS_AIR * gh);
    *new_tleaf = tair + Tdiff / 4.;
    *Cs = Ca - anleaf / gbc;                /* CO2 conc at the leaf surface */
    *dleaf = *et * press / gv;              /* VPD at the leaf surface */


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
        grashof = 1.6E8 * fabs(tleaf - tair) * leaf_width_cubed;
        gbh = 0.5 * DHEAT * pow(grashof, 0.25) / leaf_width * cmolar;
    } else {
        gbh = 0.0;
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

void update_daily_carbon_fluxes(fluxes *f, params *p, double acanopy) {

    /* umol m-2 s-1 -> gC m-2 30 min-1 */
    f->gpp_gCm2 += acanopy * UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_HLFHR;
    f->npp_gCm2 += f->gpp_gCm2 * p->cue;
    f->gpp += f->gpp_gCm2 * GRAM_C_2_TONNES_HA;
    f->npp += f->npp_gCm2 * GRAM_C_2_TONNES_HA;
    f->auto_resp += f->gpp - f->npp;

    return;
}
