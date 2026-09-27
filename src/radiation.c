
#include "radiation.h"

#define LAI_THRESH 0.001
#define RAD_THRESH 0.001
#define VIS 0
#define NIR 1

void get_diffuse_frac(canopy_wk *cw, int doy, double sw_rad) {
    /*
        For the moment, I am only going to implement Spitters, so this is a bit
        of a useless wrapper function.

    */
    spitters(cw, doy, sw_rad);

    return;
}

void spitters(canopy_wk *cw, int doy, double sw_rad) {

    /*
        Spitters algorithm to estimate the diffuse component from the measured
        irradiance.

        Eqn 20a-d.

        Parameters:
        ----------
        doy : int
            day of year
        sw_rad : double
            total incident radiation [J m-2 s-1]

        Returns:
        -------
        diffuse : double
            diffuse component of incoming radiation (returned in cw structure)

        References:
        ----------
        * Spitters, C. J. T., Toussaint, H. A. J. M. and Goudriaan, J. (1986)
          Separating the diffuse and direct component of global radiation and
          its implications for modeling canopy photosynthesis. Part I.
          Components of incoming radiation. Agricultural Forest Meteorol.,
          38:217-229.
    */
    double solar_constant, tmpr, tmpk, tmprat;

    solar_constant = 1370.0; // W m–2

    cw->direct_frac = 0.0;
    tmpr = 0.847 + cw->cos_zenith * (1.04 * cw->cos_zenith - 1.61);
    tmpk = (1.47 - tmpr) / 1.66;

    if (cw->cos_zenith > 1.0e-10 && sw_rad > 10.0) {
        tmprat = sw_rad / (solar_constant * (1.0 + 0.033 * \
                    cos(2. * M_PI * (doy-10.0) / 365.0)) * cw->cos_zenith);
    } else {
        tmprat = 0.0;
    }

    if (tmprat > 0.22) {
        cw->direct_frac = 6.4 * ( tmprat - 0.22 ) * ( tmprat - 0.22 );
    }

    if (tmprat > 0.35) {
        cw->direct_frac = MIN( 1.66 * tmprat - 0.4728, 1.0 );
    }

    if (tmprat > tmpk) {
        cw->direct_frac = MAX( 1.0 - tmpr, 0.0 );
    }

    if (cw->cos_zenith < 1.0e-2) {
        cw->direct_frac = 0.0;
    }

    cw->diffuse_frac = 1.0 - cw->direct_frac;

    return;

}

static double diffuse_extinction(params *p, double lai, double *kbx) {
    /*
        Extinction coefficient of diffuse radiation for a canopy with black
        leaves, integrating kb over three sky angles (15, 45, 75 degrees),
        eq 27 Kowalczyk et al. 2006 (CABLE). Also returns kb for those angles
        (kbx, needed for the canopy diffuse reflectance).
    */
    double gauss_w[3] = {0.308, 0.514, 0.178};
    double ang[3] = {15.0, 45.0, 75.0};
    double xphi1, xphi2, cosa, sum = 0.0;
    int    i;

    xphi1 = 0.5 - p->leaf_chi * (0.633 + 0.33 * p->leaf_chi);
    xphi2 = 0.877 * (1.0 - 2.0 * xphi1);
    for (i = 0; i < 3; i++) {
        cosa = cos(DEG2RAD(ang[i]));
        kbx[i] = (xphi1 + xphi2 * cosa) / cosa;
        sum += gauss_w[i] * exp(-kbx[i] * lai);
    }

    if (lai > LAI_THRESH) {
        return (-log(sum) / lai);
    } else {
        return (0.7);   // bare soil
    }
}

void calculate_absorbed_radiation(canopy_wk *cw, params *p, state *s,
                                  double sw_rad, double tair, double lwdown) {
    /*
        Calculate absorded irradiance of sunlit and shaded fractions of
        the canopy, the soil, and the big-leaf radiative conductances. All
        expressed on a ground-area basis. Follows CABLE (cable_albedo.F90 and
        cable_radiation.F90), which implements Wang and Leuning (1998).

        NB:  sin_beta == cos_zenith

        References:
        -----------
        * Wang and Leuning (1998) AFm, 91, 89-111. B3b and B4, the answer is
          identical de P & F
        * Kowalczyk et al. (2006) CSIRO Marine and Atmospheric Research
          paper 013 (CABLE).

        but see also:
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        * Dai et al. (2004) Journal of Climate, 17, 2281-2299.
    */
    double lai = s->lai, tk, flpwb, flwv, flws, lw_down, emissivity_air;
    double kbx[3], kd, gross, xphi1, xphi2, transb, transd, gr;
    double c1[2], rhoch[2], rhocdf[2], albsoil[2], k_dash_d[2], k_dash_b[2];
    double cexpk_dash_d[2], cexpk_dash_b[2], rho_td[2], rho_tb[2], rhocbm[2];
    double tau[2], refl[2], qsun[2] = {0.0, 0.0}, qsha[2] = {0.0, 0.0};
    double qcan_sun_lw = 0.0, qcan_sha_lw = 0.0, Ib, Id, sfact;
    double a1, a2, a3, a4, a5, a6;
    double gauss_w[3] = {0.308, 0.514, 0.178};
    int    b, vegetated_and_sunlit;

    tau[VIS] = p->leaf_tau_vis;
    tau[NIR] = p->leaf_tau_nir;
    refl[VIS] = p->leaf_refl_vis;
    refl[NIR] = p->leaf_refl_nir;
    vegetated_and_sunlit = (lai > LAI_THRESH) && (sw_rad > RAD_THRESH);

    // isothermal net radiation: leaves & soil at air temperature
    tk = tair + DEG_TO_KELVIN;
    flpwb = SIGMA * pow(tk, 4.0);           // black-body long-wave radiation
    flwv = LEAF_EMISSIVITY * flpwb;
    flws = SOIL_EMISSIVITY * flpwb;
    lw_down = downward_longwave(tair, lwdown);
    emissivity_air = lw_down / flpwb;

    // Ross-Goudriaan function is the ratio of the projected area of leaves
    // in the direction perpendicular to the direction of incident solar
    // radiation and the actual leaf area. Approximated as eqn 28,
    // Kowalcyk et al. 2006)
    xphi1 = 0.5 - p->leaf_chi * (0.633 + 0.33 * p->leaf_chi);
    xphi2 = 0.877 * (1.0 - 2.0 * xphi1);
    gross = xphi1 + xphi2 * cw->cos_zenith;

    // extinction coefficient of direct beam radiation for a canopy with black
    // leaves, eq 26 Kowalcyk et al. 2006. As CABLE this depends on the sun
    // angle only (not the beam fraction), so the sunlit leaf area doesn't
    // jump when the sky becomes overcast.
    if (lai > LAI_THRESH && cw->cos_zenith > 1.0e-6) {
        cw->kb = gross / cw->cos_zenith;
    } else {
        cw->kb = 0.5;   // bare soil
    }

    kd = diffuse_extinction(p, lai, kbx);
    if (fabs(cw->kb - kd) < RAD_THRESH) {
        cw->kb = kd + RAD_THRESH;
    }
    if (cw->cos_zenith < 1.0e-6) {
        cw->kb = 1.e5;
    }
    cw->kd = kd;

    transb = exp(-MIN(cw->kb * lai, 30.0));
    transd = lai > LAI_THRESH ? exp(-kd * lai) : 1.0;

    // soil reflectance, the soil is darker in the visible than the NIR,
    // as CABLE's albsoilsn(:,1:2)
    if (p->soil_refl <= 0.14) {
        sfact = 0.5;
    } else if (p->soil_refl <= 0.20) {
        sfact = 0.62;
    } else {
        sfact = 0.68;
    }
    albsoil[NIR] = 2.0 * p->soil_refl / (1. + sfact);
    albsoil[VIS] = sfact * albsoil[NIR];

    for (b = 0; b < 2; b++) {
        c1[b] = sqrt(1. - tau[b] - refl[b]);

        // Canopy reflection black horiz leaves
        // (eq. 6.19 in Goudriaan and van Laar, 1994):
        rhoch[b] = (1.0 - c1[b]) / (1.0 + c1[b]);

        // Canopy reflection of diffuse radiation for black leaves:
        rhocdf[b] = rhoch[b] * 2. * (gauss_w[0] * kbx[0] / (kbx[0] + kd) +
                                     gauss_w[1] * kbx[1] / (kbx[1] + kd) +
                                     gauss_w[2] * kbx[2] / (kbx[2] + kd));

        // Update extinction coefficients and fractional transmittance for
        // leaf transmittance and reflection (ie. NOT black leaves):
        // modified k diffuse(6.20)(for leaf scattering)
        k_dash_d[b] = kd * c1[b];
        cexpk_dash_d[b] = exp(-k_dash_d[b] * lai);

        // effective canopy-soil diffuse reflectance (fraction)
        if (lai > LAI_THRESH) {
            rho_td[b] = rhocdf[b] + (albsoil[b] - rhocdf[b]) *
                            (cexpk_dash_d[b] * cexpk_dash_d[b]);
        } else {
            rho_td[b] = albsoil[b];
        }

        if (vegetated_and_sunlit) {
            k_dash_b[b] = cw->kb * c1[b];
        } else {
            k_dash_b[b] = 1.e-9;
        }

        // Canopy reflection (6.21) beam:
        rhocbm[b] = 2. * cw->kb / (cw->kb + kd) * rhoch[b];

        // Canopy beam transmittance (fraction):
        cexpk_dash_b[b] = exp(-MIN(k_dash_b[b] * lai, 30.));

        // effective canopy-soil beam reflectance (fraction):
        rho_tb[b] = rhocbm[b] + (albsoil[b] - rhocbm[b]) *
                        (cexpk_dash_b[b] * cexpk_dash_b[b]);
    }

    cw->qssabs = 0.0;
    if (vegetated_and_sunlit) {
        for (b = 0; b < 2; b++) {
            // Beam and diffuse irradiance *per waveband*, shortwave is split
            // equally between the visible and NIR (as in CABLE).
            Ib = 0.5 * sw_rad * cw->direct_frac;
            Id = 0.5 * sw_rad * cw->diffuse_frac;

            // Radiation absorbed by the sunlit leaf, B3b Wang and Leuning
            // 1998: scattered diffuse, scattered beam, direct beam
            a1 = Id * (1.0 - rho_td[b]) * k_dash_d[b];
            a2 = psi_func(k_dash_d[b] + cw->kb, lai);
            a3 = Ib * (1.0 - rho_tb[b]) * k_dash_b[b];
            a4 = psi_func(k_dash_b[b] + cw->kb, lai);
            a5 = Ib * (1.0 - tau[b] - refl[b]) * cw->kb;
            a6 = psi_func(cw->kb, lai) - psi_func(2.0 * cw->kb, lai);
            qsun[b] = a1 * a2 + a3 * a4 + a5 * a6;

            // Radiation absorbed by the shaded leaf, B4  Wang and Leuning 1998
            a2 = psi_func(k_dash_d[b], lai) -
                    psi_func(k_dash_d[b] + cw->kb, lai);
            a4 = psi_func(k_dash_b[b], lai) -
                    psi_func(k_dash_b[b] + cw->kb, lai);
            qsha[b] = a1 * a2 + a3 * a4 - a5 * a6;

            // absorbed by the soil (CABLE qssabs)
            cw->qssabs += Ib * (1.0 - rho_tb[b]) * cexpk_dash_b[b] +
                          Id * (1.0 - rho_td[b]) * cexpk_dash_d[b];
        }
    } else {
        cw->qssabs = 0.5 * sw_rad * ((1.0 - albsoil[VIS]) +
                                     (1.0 - albsoil[NIR]));
    }

    if (lai > LAI_THRESH) {
        // Isothermal long-wave absorbed by the sunlit & shaded leaves (CABLE
        // qcan(:,:,3)): exchange with the soil, the sky and the other leaves
        qcan_sun_lw = (flws - flwv) * kd * (transd - transb) /
                        (cw->kb - kd) +
                      (emissivity_air - LEAF_EMISSIVITY) * kd * flpwb *
                        (1.0 - transd * transb) / (cw->kb + kd);
        qcan_sha_lw = (1.0 - transd) * (flws + lw_down - 2.0 * flwv) -
                        qcan_sun_lw;

        // Radiative conductance of the big leaves (mol m-2 s-1), CABLE
        // gradis. A leaf deep in the canopy mostly sees other leaves at
        // the same temperature, so this is much less than 2 x the single
        // leaf value x LAI.
        gr = 4.0 * LEAF_EMISSIVITY * SIGMA * tk * tk * tk / (CP * MASS_AIR);
        cw->gradis[SUNLIT] = gr * kd * ((1.0 - transb * transd) /
                                        (cw->kb + kd) +
                                        (transd - transb) / (cw->kb - kd));
        cw->gradis[SHADED] = 2.0 * gr * (1.0 - transd) - cw->gradis[SUNLIT];
        cw->gradis[SUNLIT] = MAX(1.0e-3, cw->gradis[SUNLIT]);
        cw->gradis[SHADED] = MAX(1.0e-3, cw->gradis[SHADED]);
    } else {
        cw->gradis[SUNLIT] = 1.0e-3;
        cw->gradis[SHADED] = 1.0e-3;
    }

    // soil net radiation: SW reaching the soil, LW from the sky through the
    // gaps and from the canopy, minus soil emission
    cw->rnet_soil = cw->qssabs + transd * lw_down + (1.0 - transd) * flwv -
                    flws;

    cw->apar_leaf[SUNLIT] = qsun[VIS] * J_2_UMOL;
    cw->apar_leaf[SHADED] = qsha[VIS] * J_2_UMOL;

    // Total energy absorbed by canopy, summing VIS, NIR and LW components, to
    // leave us with the indivual leaf components.
    cw->rnet_leaf[SUNLIT] = qsun[VIS] + qsun[NIR] + qcan_sun_lw;
    cw->rnet_leaf[SHADED] = qsha[VIS] + qsha[NIR] + qcan_sha_lw;

    if (vegetated_and_sunlit) {
        /* Calculate sunlit &shdaded LAI of the canopy - de P * F eqn 18*/
        cw->lai_leaf[SUNLIT] = (1.0 - transb) / cw->kb;
        cw->lai_leaf[SHADED] = lai - cw->lai_leaf[SUNLIT];
    } else {
        cw->lai_leaf[SUNLIT] = 0.0;
        cw->lai_leaf[SHADED] = 0.0;
    }

    return;
}

void calculate_soil_net_radiation_night(canopy_wk *cw, params *p, state *s,
                                        double tair, double lwdown) {
    /* soil net (long-wave) radiation when the sun is down, as above */
    double kbx[3], transd, flpwb;

    flpwb = SIGMA * pow(tair + DEG_TO_KELVIN, 4.0);
    transd = s->lai > LAI_THRESH ? exp(-diffuse_extinction(p, s->lai, kbx) *
                                       s->lai) : 1.0;
    cw->qssabs = 0.0;
    cw->rnet_soil = transd * downward_longwave(tair, lwdown) +
                    (1.0 - transd) * LEAF_EMISSIVITY * flpwb -
                    SOIL_EMISSIVITY * flpwb;

    return;
}

double downward_longwave(double tair, double lwdown) {
    /*
        Downward long-wave (W m-2): the forcing if it has it, otherwise the
        clear sky estimate from air temperature (K), Swinbank, W. C. (1963)
        Q. J. R. Meteorol. Soc., 89, 339–348.
    */
    double tk = tair + DEG_TO_KELVIN;

    if (lwdown > 0.0) {
        return (lwdown);
    }
    return (0.0000094 * SIGMA * pow(tk, 6.0));
}

double psi_func(double z, double lai) {
    /*
        B5 function from Wang and Leuning which integrates property passed via
        arg list over the canopy space

        References:
        -----------
        * Wang and Leuning (1998) AFm, 91, 89-111. Page 106

    */
    return ( (1.0 - exp(-MIN(z * lai, 30.0))) / z );
}


void calculate_solar_geometry(canopy_wk *cw, params *p, double doy,
                              double hod) {

    /*
        The solar zenith angle is the angle between the zenith and the centre
        of the sun's disc. The solar elevation angle is the altitude of the
        sun, the angle between the horizon and the centre of the sun's disc.
        Since these two angles are complementary, the cosine of either one of
        them equals the sine of the other, i.e. cos theta = sin beta. I will
        use cos_zen throughout code for simplicity.

        Arguments:
        ----------
        params : p
            params structure
        doy : double
            day of year
        hod : double:
            hour of the day [0.5 to 24]
        cos_zen : double
            cosine of the zenith angle of the sun in degrees (returned)
        elevation : double
            solar elevation (degrees) (returned)

        References:
        -----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.

    */

    double sindec, zenith_angle;

    /* need to convert 30 min data, 0-47 to 0-23.5 */
    hod /= 2.0;

    // sine of maximum declination
    sindec = -sin(23.45 * M_PI / 180.) * \
                cos(2. * M_PI * (doy + 10.0) / 365.0);

    cw->cos_zenith = MAX(sin(M_PI / 180. * p->latitude) * sindec + \
                         cos(M_PI / 180. * p->latitude) * \
                         sqrt(1. - sindec * sindec) * \
                         cos(M_PI * (hod - 12.0) / 12.0), 1e-8);

    zenith_angle = RAD2DEG(acos(cw->cos_zenith));
    cw->elevation = 90.0 - zenith_angle;

    return;
}
