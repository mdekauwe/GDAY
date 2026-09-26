
#include "radiation.h"

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

void calculate_absorbed_radiation(canopy_wk *cw, params *p, state *s,
                                  double sw_rad, double tair) {
    /*
        Calculate absorded irradiance of sunlit and shaded fractions of
        the canopy. The total irradiance absorbed by the canopy and the
        sunlit/shaded components are all expressed on a ground-area basis!

        NB:  sin_beta == cos_zenith

        References:
        -----------
        * Wang and Leuning (1998) AFm, 91, 89-111. B3b and B4, the answer is
          identical de P & F

        but see also:
        * De Pury & Farquhar (1997) PCE, 20, 537-557.
        * Dai et al. (2004) Journal of Climate, 17, 2281-2299.
    */

    double lw_down = 0.0, flpwb = 0.0, emissivity_air = 0.0;
    double cos3_15 = 0.0, cos3_45 = 0.0, cos3_75 = 0.0;
    double tk = 0.0, txx1 = 0.0, txx2 = 0.0;
    double txx3 = 0.0, kbx1 = 0.0, kbx2 = 0.0, kbx3 = 0.0;
    double c1_1 = 0.0, c1_2 = 0.0, rhoch_1 = 0.0, rhoch_2 = 0.0;
    double rhocdf_vis= 0.0, rhocdf_nir= 0.0, sfact= 0.0;
    double albsoil_vis = 0.0, albsoil_nir = 0.0;
    double k_dash_d_vis = 0.0, k_dash_d_nir = 0.0, cexpk_dash_d_vis = 0.0;
    double cexpk_dash_d_nir = 0.0, k_dash_b_nir = 0.0, cexpk_dash_b_nir = 0.0;
    double rho_td_vis = 0.0, rho_td_nir = 0.0, k_dash_b_vis = 0.0;
    double rhocbm_vis = 0.0, rhocbm_nir = 0.0, cexpk_dash_b_vis = 0.0;
    double rho_tb_vis = 0.0, rho_tb_nir = 0.0;
    double a1_vis = 0.0, a2_vis = 0.0, a3_vis = 0.0, a4_vis = 0.0;
    double a5_vis = 0.0, a6_vis = 0.0;
    double a1_nir = 0.0, a2_nir = 0.0, a3_nir = 0.0, a4_nir = 0.0;
    double a5_nir = 0.0, a6_nir = 0.0;
    double qcan_sha_vis = 0.0, qcan_sun_vis = 0.0;
    double qcan_sha_nir = 0.0, qcan_sun_nir = 0.0;
    double qcan_sha_lw = 0.0, qcan_sun_lw = 0.0, Ib = 0.0, Id = 0.0;
    double xphi1 = 0.0, xphi2 = 0.0, Gross = 0.0, kd = 0.0;
    double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4 = 0.0, a5 = 0.0, a6 = 0.0;
    // leaf emissivity (-), Table 3, Wang and Leuning, 1998
    double emissivity_leaf = 0.96;

    // soil emissivity (-), Table 3, Wang and Leuning, 1998
    double emissivity_soil = 0.94;

    double LAI_THRESH = 0.001;
    double RAD_THRESH = 0.001;

    // Gaussian integ. weights
    double gauss_w1 = 0.308;
    double gauss_w2 = 0.514;
    double gauss_w3 = 0.178;

    // leaf transmissivity [-] (VIS: 0.07 - 0.15)
    // ENF: 0.05; EBF: 0.05; DBF: 0.05; C3G: 0.070
    double tau_vis = 0.1;
    double tau_nir = 0.3;

    // leaf reflectance [-] (VIS:0.07 - 0.15)
    // ENF: 0.062;EBF: 0.076;DBF: 0.092; C3G: 0.11
    double refl_vis = 0.1;
    double refl_nir = 0.3;

    // Table 3, Wang and Leuning, 1998
    double soil_reflectance = 0.1; //(same as MAESTRA)

    //empirical param related to the leaf angle dist (= 0 for spherical LAD)
    double chi = 9.99999978E-03;

    // surface temperaute - just using air temp
    tk = tair + DEG_TO_KELVIN;
    

    // Estimate LWdown based on an emprical function of air temperature (K)
    //following Swinbank, W. C. (1963): Long-wave radiation from clear skies,
    // Q. J. R. Meteorol. Soc., 89, 339–348, doi:10.1002/qj.49708938105.
    lw_down = 0.0000094 * SIGMA * pow(tk, 6.0);

    // black-body long-wave radiation
    flpwb = SIGMA * pow(tk, 4.0);

    // air emissivity
    emissivity_air = lw_down / flpwb;

    // cos(15 45 75 degrees)
    cos3_15 = cos(DEG2RAD(15.0));
    cos3_45 = cos(DEG2RAD(45.0));
    cos3_75 = cos(DEG2RAD(75.0));

    // leaf angle parmameter 1
    xphi1 = 0.5 - chi * (0.633 + 0.33 * chi);

    // leaf angle parmameter 2
    xphi2 = 0.877 * (1.0 - 2.0 * xphi1);

    // Ross-Goudriaan function is the ratio of the projected area of leaves
    // in the direction perpendicular to the direction of incident solar
    // radiation and the actual leaf area. Approximated as eqn 28,
    // Kowalcyk et al. 2006)
    Gross = xphi1 + xphi2 * cw->cos_zenith;

    // extinction coefficient of direct beam radiation for a canopy with black
    // leaves, eq 26 Kowalcyk et al. 2006
    if ( (s->lai > LAI_THRESH) & (cw->direct_frac > RAD_THRESH) ) {   // vegetated
        cw->kb = Gross / cw->cos_zenith;
    } else {   // i.e. bare soil
        cw->kb = 0.5;
    }

    // extinction coefficient of diffuse radiation for a canopy with black
    // leaves, eq 27 Kowalcyk et al. 2006
    if (s->lai > LAI_THRESH) {  // vegetated

        // Approximate integration of kb
        kbx1 = (xphi1 + xphi2 * cos3_15) / cos3_15;
        kbx2 = (xphi1 + xphi2 * cos3_45) / cos3_45;
        kbx3 = (xphi1 + xphi2 * cos3_75) / cos3_75;

        txx1 = gauss_w1 * exp(-kbx1 * s->lai);
        txx2 = gauss_w2 * exp(-kbx2 * s->lai);
        txx3 = gauss_w3 * exp(-kbx3 * s->lai);

        kd = -log(txx1 + txx2 + txx3) / s->lai;
    } else {   // i.e. bare soil
        kd = 0.7;
    }

    if (fabs(cw->kb - kd) < RAD_THRESH) {
        cw->kb = kd + RAD_THRESH;
    }

    if (cw->direct_frac < RAD_THRESH) {
        cw->kb = 1.e5;
    }

    c1_1 = sqrt(1. - tau_vis - refl_vis);
    c1_2 = sqrt(1. - tau_nir - refl_nir);

    // Canopy reflection black horiz leaves
    // (eq. 6.19 in Goudriaan and van Laar, 1994):
    rhoch_1 = (1.0 - c1_1) / (1.0 + c1_1);
    rhoch_2 = (1.0 - c1_2) / (1.0 + c1_2);

    // Canopy reflection of diffuse radiation for black leaves:
    rhocdf_vis = rhoch_1 * 2. * \
                    (gauss_w1 * kbx1 / (kbx1 + kd) + \
                     gauss_w2 * kbx2 / (kbx2 + kd) + \
                     gauss_w3 * kbx3 / (kbx3 + kd));

    rhocdf_nir = rhoch_2 * 2. * \
                    (gauss_w1 * kbx1 / (kbx1 + kd) + \
                     gauss_w2 * kbx2 / (kbx2 + kd) + \
                     gauss_w3 * kbx3 / (kbx3 + kd));

    // Calculate albedo
    if (soil_reflectance <= 0.14) {
        sfact = 0.5;
    } else if ( ( soil_reflectance > 0.14) & (soil_reflectance <= 0.20) ) {
        sfact = 0.62;
    } else {
        sfact = 0.68;
    }

    // soil + snow reflectance (ignoring snow), as CABLE's albsoilsn(:,1:2),
    // i.e. the soil is darker in the visible than the NIR
    albsoil_nir = 2.0 * soil_reflectance / (1. + sfact);
    albsoil_vis = sfact * albsoil_nir;

    // Update extinction coefficients and fractional transmittance for
    // leaf transmittance and reflection (ie. NOT black leaves):
    // modified k diffuse(6.20)(for leaf scattering)
    k_dash_d_vis = kd * c1_1;
    k_dash_d_nir = kd * c1_2;

    // Define canopy diffuse transmittance (fraction):
    cexpk_dash_d_vis = exp(-k_dash_d_vis * s->lai);
    cexpk_dash_d_nir = exp(-k_dash_d_nir * s->lai);

    // Calculate effective canopy-soiil diffuse reflectance (fraction)
    if (s->lai > 0.001) {
        rho_td_vis = rhocdf_vis + (albsoil_vis - rhocdf_vis) * \
                        (cexpk_dash_d_vis * cexpk_dash_d_vis);
        rho_td_nir = rhocdf_nir + (albsoil_nir - rhocdf_nir) * \
                        (cexpk_dash_d_nir * cexpk_dash_d_nir);
    } else {
        rho_td_vis = albsoil_vis;
        rho_td_nir = albsoil_nir;
    }

    // where vegetated and sunlit
    if ( (s->lai > LAI_THRESH) & (sw_rad > RAD_THRESH) ) {
        k_dash_b_vis = cw->kb * c1_1;
        k_dash_b_nir = cw->kb * c1_2;
    } else {
        k_dash_b_vis = 1.e-9;
        k_dash_b_nir = 1.e-9;
    }

    // Canopy reflection (6.21) beam:
    rhocbm_vis = 2. * cw->kb / (cw->kb + kd) * rhoch_1;
    rhocbm_nir = 2. * cw->kb / (cw->kb + kd) * rhoch_2;

    // Canopy beam transmittance (fraction):
    cexpk_dash_b_vis = exp(-MIN(k_dash_b_vis * s->lai, 30.));
    cexpk_dash_b_nir = exp(-MIN(k_dash_b_nir * s->lai, 30.));

    // Calculate effective canopy-soil beam reflectance (fraction):
    rho_tb_vis = rhocbm_vis + (albsoil_vis - rhocbm_vis) * \
                    (cexpk_dash_b_vis * cexpk_dash_b_vis);
    rho_tb_nir = rhocbm_nir + (albsoil_nir - rhocbm_nir) * \
                    (cexpk_dash_b_nir * cexpk_dash_b_nir);

    if ( (s->lai > LAI_THRESH) & (sw_rad > RAD_THRESH) ) {

        // Beam and diffuse irradiance *per waveband*, shortwave is split
        // equally between the visible and NIR (as in CABLE). Using the full
        // shortwave for each band absorbed ~2x the incident PAR.
        Ib = 0.5 * sw_rad * cw->direct_frac;
        Id = 0.5 * sw_rad * cw->diffuse_frac;

        a1_vis = Id * (1.0 - rho_td_vis) * k_dash_d_vis;
        a1_nir = Id * (1.0 - rho_td_nir) * k_dash_d_nir;

        a2_vis = psi_func(k_dash_d_vis + cw->kb, s->lai);
        a2_nir = psi_func(k_dash_d_nir + cw->kb, s->lai);

        a3_vis = Ib * (1.0 - rho_tb_vis) * k_dash_b_vis;
        a3_nir = Ib * (1.0 - rho_tb_nir) * k_dash_b_nir;

        a4_vis = psi_func(k_dash_b_vis + cw->kb, s->lai);
        a4_nir = psi_func(k_dash_b_nir + cw->kb, s->lai);

        a5_vis = Ib * (1.0 - tau_vis - refl_vis) * cw->kb;
        a5_nir = Ib * (1.0 - tau_nir - refl_nir) * cw->kb;

        a6_vis = psi_func(cw->kb, s->lai) - psi_func(2.0 * cw->kb, s->lai);
        a6_nir = psi_func(cw->kb, s->lai) - psi_func(2.0 * cw->kb, s->lai);

        qcan_sun_vis = (a1_vis * a2_vis + a3_vis * a4_vis + a5_vis * a6_vis);
        qcan_sun_nir = (a1_nir * a2_nir + a3_nir * a4_nir + a5_nir * a6_nir);

        // Radiation absorbed by the shaded leaf, B4  Wang and Leuning 1998
        a2_vis = psi_func(k_dash_d_vis, s->lai) - \
                    psi_func(k_dash_d_vis + cw->kb, s->lai);
        a2_nir = psi_func(k_dash_d_nir, s->lai) - \
                    psi_func(k_dash_d_nir + cw->kb, s->lai);

        a4_vis = psi_func(k_dash_b_vis, s->lai) - \
                    psi_func(k_dash_b_vis + cw->kb, s->lai);
        a4_nir = psi_func(k_dash_b_nir, s->lai) - \
                    psi_func(k_dash_b_nir + cw->kb, s->lai);

        qcan_sha_vis = (a1_vis * a2_vis + a3_vis * a4_vis - a5_vis * a6_vis);
        qcan_sha_nir = (a1_nir * a2_nir + a3_nir * a4_nir - a5_nir * a6_nir);

    }

    // Longwave radiation absorbed by sunlit leaves under isothermal conditions
    // B18 Wang and Leuning 1998
    a1 = -kd * SIGMA * pow(tk, 4);
    a2 = emissivity_leaf * (1.0 - emissivity_air);
    a3 = psi_func(cw->kb + kd, s->lai);
    a4 = 1.0 - emissivity_soil;
    a5 = (emissivity_leaf - emissivity_air);
    a6 = psi_func(2.0 * kd, s->lai) * psi_func(cw->kb - kd, s->lai);
    qcan_sun_lw = a1 * (a2 * a3 + a4 * a5 * a6);

    // Longwave radiation absorbed by shaded leaves under isothermal conditions
    // B19 Wang and Leuning 1998
    a3 = psi_func(kd, s->lai);
    a6 = exp(-kd * s->lai) * a3;
    qcan_sha_lw = a1 * (a2 * a3 - a4 * a5 * a6) - qcan_sun_lw;

    cw->apar_leaf[SUNLIT] = qcan_sun_vis * J_2_UMOL;
    cw->apar_leaf[SHADED] = qcan_sha_vis * J_2_UMOL;

    // Total energy absorbed by canopy, summing VIS, NIR and LW components, to
    // leave us with the indivual leaf components.
    cw->rnet_leaf[SUNLIT] = qcan_sun_vis + qcan_sun_nir + qcan_sun_lw;
    cw->rnet_leaf[SHADED] = qcan_sha_vis + qcan_sha_nir + qcan_sha_lw;

    // where vegetated and sunlit
    if ( (s->lai > LAI_THRESH) & (sw_rad > RAD_THRESH) ) {

        /* Calculate sunlit &shdaded LAI of the canopy - de P * F eqn 18*/
        cw->lai_leaf[SUNLIT] = (1.0 - exp(-cw->kb * s->lai)) / cw->kb;
        cw->lai_leaf[SHADED] = s->lai - cw->lai_leaf[SUNLIT];
    } else {
        cw->lai_leaf[SUNLIT] = 0.0;
        cw->lai_leaf[SHADED] = 0.0;
    }


    return;
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
