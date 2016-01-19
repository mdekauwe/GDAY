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

    double Cs, dleaf, tleaf, new_tleaf, et, leafn, fc, ncontent, diffuse_frac,
           zenith_angle, elevation, anleaf, gsc, apar, anir, athermal
           direct_apar, diffuse_apar, arg1, arg2;
    int    hod, iter = 0, itermax = 100, leaf;
    long   offset;
    double gpp_gCm2_30_min, SEC_2_30min = 1800.;

    /* initialise values of leaf temp, leaf surface CO2 and VPD from air space*/
    tleaf = m->tair[offset];            /* Leaf temperature */
    dleaf = m->vpd[offset] * KPA_2_PA;  /* VPD at the leaf suface */
    Cs = m->co2[offset];                /* CO2 conc. at the leaf suface */


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

    /* zero daily fluxes */
    f->gpp_gCm2 = 0.0;
    f->npp_gCm2 = 0.0;
    f->gpp = 0.0;
    f->npp = 0.0;
    f->auto_resp = 0.0;
    f->apar = 0.0;

    gpp_gCm2_30_min = 0.0;

    for (hod = 0; hod < c->num_hlf_hrs; hod++) {
        offset = project_day * c->num_days + hod;

        zenith_angle = calculate_zenith_angle(p, m->doy[project_day], hod);
        elevation = 90.0 - zenith_angle;

        /* calculates diffuse frac from half-hourly incident radiation */
        diffuse_frac = get_diffuse_frac(m->doy[offset], zenith_angle,
                                        m->par[offset]);

        /*printf("%d %lf\n", hod, elevation); */

        /* Is the sun up? If so calculate photosynthesis */
        if (elevation > 0.0) {

            while (TRUE) {

                calculate_absorbed_radiation(p, diffuse_frac, &direct_apar,
                                             &diffuse_apar);

                /* sunlit, shaded loop */
                for (leaf = 0; leaf <= 1; leaf++) {





                    if (c->ps_pathway == C3) {
                        photosynthesis_C3(c, m, p, s, offset, ncontent, tleaf,
                                          Cs, dleaf, &gsc, &anleaf);
                    } else {
                        /* Nothing implemented */
                        fprintf(stderr, "C4 photosynthesis not implemented\n");
                        exit(EXIT_FAILURE);
                    }

                    if (anleaf > 0.0) {
                        /* Calculate new Cs, dleaf, Tleaf */
                        solve_leaf_energy_balance(f, m, p, s, offset, tleaf,
                                                  gsc, anleaf, apar, anir,
                                                  athermal, &Cs, &dleaf,
                                                  &new_tleaf, &et);
                    } else {
                        et = 0.0;
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
            }


            /* umol m-2 s-1 -> gC m-2 30 min-1 */
            gpp_gCm2_30_min += anleaf * UMOL_TO_MOL * MOL_C_TO_GRAMS_C * SEC_2_30min;
            printf("* %lf %lf\n", m->par[offset], anleaf);
            exit(1);


            /* Calculate plant respiration */
            if (c->respiration_model == FIXED) {
                /* Plant respiration assuming carbon-use efficiency. */
                f->auto_resp += f->gpp * p->cue;
            } else if(c->respiration_model == TEMPERATURE) {
                fprintf(stderr, "Not implemented yet\n");
                exit(EXIT_FAILURE);
            } else if (c->respiration_model == BIOMASS) {
                fprintf(stderr, "Not implemented yet\n");
                exit(EXIT_FAILURE);
            }

            /* Calculate NPP */
            f->npp_gCm2 = f->gpp_gCm2 * p->cue;
            f->npp = f->npp_gCm2 * GRAM_C_2_TONNES_HA;


        } else {
            /* set time slot photosynthesis/respiration to be zero, but we
               still need to calc the full water balance */

        }

    }
    exit(1);

    return;

}

void calculate_absorbed_radiation(params *p, double diffuse_frac,
                                  double *direct_apar, double *diffuse_apar) {
    /*
        Calculate beam, diffuse and scattered radiation
    */

    /* Calculate beam radiation absorbed by sunlit leaf area. */
    arg1 = m->par[offset] * (1.0 - diffuse_frac);
    arg2 = cos(zenith_angle) * p->leaf_abs;
    *direct_apar = arg1 / arg2;

    /*  Calculate diffuse radiation absorbed directly. */
    *diffuse_apar = m->par[offset] * diffuse_frac;



    return;
}

void solve_leaf_energy_balance(fluxes *f, met *m, params *p, state *s,
                               long offset, double tleaf, double gsc,
                               double anleaf, double apar, double anir,
                               double athermal, double *Cs, double *dleaf,
                               double *new_tleaf, double *et) {
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
    double par = m->par[offset];
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

    Tk = tair + DEG_TO_KELVIN;

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
    sw_rad = par * PAR_2_SW; /* W m-2 */

    /* isothermal net LW radiaiton at top of canopy, assuming emissivity of
       the canopy is 1 */
    net_lw_rad = (1.0 - emissivity_atm) * SIGMA * pow(Tk, 4.0);
    rnet = p->leaf_abs * sw_rad - net_lw_rad * kd * exp(-kd * s->lai);

    /*rnet = (apar / UMOLPERJ) + anir + athermal;*/

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

    printf("** %lf %lf %lf %lf %lf %lf\n", *et, tleaf, *new_tleaf, *dleaf, *Cs, rnet);


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
