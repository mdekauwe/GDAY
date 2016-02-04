
#include "radiation.h"

double get_diffuse_frac(int doy, double cos_zenith, double par) {
    /*
        For the moment, I am only going to implement Spitters, so this is a bit
        of a useless wrapper function.

    */

    return spitters(doy, cos_zenith, par);
}

double spitters(int doy, double cos_zenith, double par) {

    /*

    Spitters algorithm to estimate the diffuse component from the measured
    irradiance.

    Parameters:
    ----------
    doy : int
        day of year
    cos_zenith : double
        sun zenith angle [radians]
    par : double
        total par measured [umol m-2 s-1]

    Returns:
    -------
    diffuse : double
        diffuse component of incoming radiation

    References:
    ----------
    * Spitters, C. J. T., Toussaint, H. A. J. M. and Goudriaan, J. (1986)
      Separating the diffuse and direct component of global radiation and its
      implications for modeling canopy photosynthesis. Part I. Components of
      incoming radiation. Agricultural Forest Meteorol., 38:217-229.
    */

    double sw_rad, sin_beta, So, tau, R, K, diffuse_frac;
    double SEC_TO_HFHR = 60.0 * 30.0;
    double solar_constant;

    /* SW_down [W/m2] = [J m-2 s-1] */
    sw_rad = par * PAR_2_SW;

    /* sine of the elevation of the sun above the horizon */
    sin_beta = DEG2RAD(90.0) - cos_zenith;
    So = calc_extra_terrestrial_irradiance(doy, sin_beta);

    /* atmospheric transmisivity */
    tau = estimate_clearness(sw_rad, So);

    /* For zenith angles > 80 degrees, diffuse_frac = 1.0 */
    if (sin_beta > 0.17) {
        /* the ratio between diffuse and total Solar irradiance (R), eqn 20 */
        /*R = 0.847 - 1.61 * cos_zenith + 1.04 * (cos_zenith * cos_zenith); */
        R = 0.847 - 1.61 * sin_beta + 1.04 * (sin_beta * sin_beta);
        K = (1.47 - R) / 1.66;
        if (tau <= 0.22) {
            diffuse_frac = 1.0;
        } else if (tau <= 0.35) {
            diffuse_frac = 1.0 - 6.4 * (tau - 0.22) * (tau - 0.22);
        } else if (tau <= K) {
            diffuse_frac = 1.47 - 1.66 * tau;
        } else {
            diffuse_frac = R;
        }
    } else {
        diffuse_frac = 1.0;
    }

    if (diffuse_frac <= 0.0) {
        diffuse_frac = 0.0;
    } else if (diffuse_frac >= 1.0) {
        diffuse_frac = 1.0;
    }

    return (diffuse_frac);

}

void calculate_absorbed_radiation(params *p, state *s, double par,
                                  double diffuse_frac, double elevation,
                                  double cos_zenith, double *apar,
                                  double *sunlit_lai, double *shaded_lai) {
    /*
        Calculate absorded irradiance of sunlit and shaded fractions of
        the canopy

        References:
        -----------
        * De Pury & Farquhar (1997) PCE, 20, 537-557.

        but see also:
        * Wang and Leuning (1998) AFm, 91, 89-111.
        * Dai et al. (2004) Journal of Climate, 17, 2281-2299.
    */

    int    i;
    double czen, integral, kb, kd, phi_1, phi_2, Gross, psi, Ib, Id, Is, Ic,
           k_dash_b, k_dash_d, scattered, shaded, beam;
    double direct_frac = 1.0 - diffuse_frac;
    double lai = s->lai;
    double lad = p->lad;

    /* canopy reflection coeffcient for diffuse PAR; de Pury & Farquhar, 1997 */
    double rho_cd = 0.036;

    /* canopy reflection coeffcient for direct PAR; de Pury & Farquhar, 1997 */
    double rho_cb = 0.029;

    /* leaf scattering coefficient of PAR; de Pury & Farquhar, 1997 */
    double omega_PAR = 0.15;

    /* direct PAR extinction coefficent - Dai et al 2004, eqn 2. */
    /*phi_1 = 0.5 - (0.633 * lad) - (0.33 * lad * lad);
    phi_2 = 0.877 * (1.0 - 2.0 * phi_1);
    Gross = phi_1 + (phi_2 * cos_zenith);
    *kb = Gross / cos_zenith;
    */

    /* beam radiation extinction coefficent of canopy - de P & Far '97, Tab 3 */
    kb = 0.5 / sin(DEG2RAD(elevation));

    /* beam & scattered PAR extinction coefficent - de P & Farq '97, Table 3*/
    k_dash_b = 0.46 / sin(DEG2RAD(elevation));

    /* diffuse & scattered PAR extinction coeff - de P & Farq '97, Table 3 */
    k_dash_d = 0.718;

    /* Diffuse beam irradiance - de Pury & Farquhar (1997), eqn 20c */
    Id = par * diffuse_frac;
    shaded = ( Id * (1.0 - rho_cd) *
                    (1.0 - exp(-(k_dash_d + kb) * lai)) *
                    (k_dash_d / (k_dash_d + kb)) );

    /* Direct beam irradiance - de Pury & Farquhar (1997), eqn 20b */
    Ib = par * direct_frac;
    beam = Ib * (1.0 - omega_PAR) * (1.0 - exp(-kb * lai));

    /* scattered-beam irradiance - de Pury & Farquhar (1997), eqn 20d */
    scattered = Ib * ( (1.0 - rho_cb) * (1.0 - exp(-(k_dash_b + kb) * lai)) *
                        k_dash_b / (k_dash_b + kb) - (1.0 - omega_PAR) *
                        (1.0 - exp(-2.0 * kb * lai)) / 2.0 );

    /* Irradiance absorbed by the canopy - de Pury & Farquhar (1997), eqn 13 */
    Ic = ( (1.0 - rho_cb) * Ib * (1.0 - exp(-k_dash_b * lai)) +
           (1.0 - rho_cd) * Id * (1.0 - exp(-k_dash_d * lai)) );

    /*
        Irradiance absorbed by the sunlit fraction of the canopy is the sum of
        direct-beam, diffuse and scattered-beam components
    */
    *(apar+SUNLIT) = beam + scattered + shaded;

    /*
        Irradiance absorbed by the shaded leaf area of the canopy is the
        integral of absorbed irradiance in the shade (Eqn A7) and the shdaded
        leaf area fraction. Or simply, the difference between the total
        irradiacne absorbed by the canopy and the irradiance absorbed by the
        sunlit leaf area
    */
    *(apar+SHADED) = Ic - *(apar+SUNLIT);

    /*
        Scale leaf fluxes to the canopy
        - Fractional area decreases exponentialy with LAI from the top
          of the canopy. Integrating sunlit/shaded fraction over the
          canopy to calculate leaf area index of sunlit/shaded fractions
          of the canopy. De Pury & Farquhar 1997, eqn 18.
    */
    *sunlit_lai = (1.0 - exp(-kb * s->lai)) / kb;
    *shaded_lai = s->lai - *sunlit_lai;

    
    /*printf("%lf %lf %lf\n", par, *(apar+SUNLIT), *(apar+SHADED));
    exit(1);*/
    return;
}

void calculate_zenith_angle(params *p, double doy, double hod,
                            double *cos_zen, double *elevation) {

    /*
    Estimate the sun zenith angle (Degrees)

    Arguments:
    ----------
    hour -- [0.5 to 24]
    lat -- latitude in degrees
    lon -- longitude in degrees
    doy -- day of year

    Returns:
    --------
    zen_angle : float
        zenith angle of the sun in degrees
    cos_zen : float
        cosine of zenith angle

    References:
    -----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.

    */

    double dec, et, lc, t0, h, gamma, zenith_angle, hour_angle, rlat;

    /* need to convert 30 min data, 0-47 to 0-23.5 */
    hod /= 2.0;

    gamma = day_angle(doy);
    dec = calculate_solar_declination(doy, gamma);
    et = calculate_eqn_of_time(gamma);
    t0 = calculate_solar_noon(et, p->longitude);
    h = calculate_hour_angle(hod, t0);
    rlat = DEG2RAD(p->latitude);

    /* solar elevation angle (degrees) A13 - De Pury & Farquhar  */
    *elevation = (RAD2DEG(asin(sin(rlat) * sin(dec) + cos(rlat) *
                          cos(dec) * cos(h))));
    zenith_angle = 90.0 - *elevation;
    *cos_zen = cos(DEG2RAD(zenith_angle));
    if (*cos_zen > 1.0)
        *cos_zen = 1.0;
    else if (*cos_zen < 0.0)
        *cos_zen = 0.0;

    return;
}

double calculate_solar_noon(double et, double longitude) {
    /* Calculation solar noon - De Pury & Farquhar, '97: eqn A16

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.

    Returns:
    ---------
    t0 - solar noon (hours).
    */
    double t0, Ls;

    /* all international standard meridians are multiples of 15deg east/west of
       greenwich */
    Ls = round_to_value(longitude, 15.);
    t0 = 12.0 + (4.0 * (Ls - longitude) - et) / 60.0;

    return (t0);
}

double calculate_hour_angle(double t, double t0) {
    /* Calculation solar noon - De Pury & Farquhar, '97: eqn A15

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.

    Returns:
    ---------
    h - hour angle (radians).
    */
    return (M_PI * (t - t0) / 12.0);

}

double day_angle(int doy) {
    /* Calculation of day angle - De Pury & Farquhar, '97: eqn A18

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.

    Returns:
    ---------
    gamma - day angle in radians.
    */

    return (2.0 * M_PI * (doy - 1.0) / 365.0);
}

double calculate_solar_declination(int doy, double gamma) {
    /*
    Solar Declination Angle is a function of day of year and is indepenent
    of location, varying between 23deg45' to -23deg45'

    Arguments:
    ----------
    doy : int
        day of year, 1=jan 1
    gamma : double
        fractional year (radians)

    Returns:
    --------
    dec: float
        Solar Declination Angle [radians]

    Reference:
    ----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    */
    double sindec, decl;

    /* declination (radians) */
    /*decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - \
           0.006758 * cos(2.0 * gamma) + 0.000907 * sin(2.0 * gamma) -\
           0.002697 * cos(3.0 * gamma) + 0.00148 * sin(3.0 * gamma);*/


    /* (radians) A14 - De Pury & Farquhar  */
    decl = -23.4 * (M_PI / 180.) * cos(2.0 * M_PI * (doy + 10) / 365);

    return (decl);

}

double calculate_eqn_of_time(double gamma) {
    /* Equation of time - correction for the difference btw solar time
    and the clock time.

    Arguments:
    ----------
    doy : int
        day of year
    gamma : double
        fractional year (radians)

    References:
    -----------
    * De Pury & Farquhar (1997) PCE, 20, 537-557.
    * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
      biophysics. Pg 169.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    * Hughes, David W.; Yallop, B. D.; Hohenkerk, C. Y. (1989),
      "The Equation of Time", Monthly Notices of the Royal Astronomical
      Society 238: 1529–1535
    */
    double et, f, A;

    /* radians */
    /*et = 0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -\
         0.014615 * cos(2.0 * gamma) - 0.04089 * sin(2.0 * gamma);*/

    /* radians to minutes */
    /*et *= 229.18; */

    /* radians to hours */
    /*et *= 24.0 / (2.0 * M_PI);*/

    /* minutes - de Pury and Farquhar, 1997 - A17 */
    et = (0.017 + 0.4281 * cos(gamma) - 7.351 * sin(gamma) - 3.349 *
          cos(2.0 * gamma) - 9.731  * sin(gamma));

    return (et);
}



double calc_extra_terrestrial_irradiance(double doy, double sin_beta) {
    /* Solar radiation incident outside the earth's atmosphere, e.g.
    extra-terrestrial radiation. The value varies a little with the earths
    orbit.

    Using formula from Spitters not Leuning!

    Arguments:
    ----------
    doy : double
        day of year
    sin_beta : double
        sine of the elevation of the sun above the horizon (radians)

    Returns:
    --------
    So : float
        solar radiation normal to the sun's bean outside the Earth's atmosphere
        (W m-2)

    Reference:
    ----------
    * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    * Spitters et al. (1986) AFM, 38, 217-229.
    */
    double Sc = 1370.0; /* W m-2 */
    double orbit_correction;

    /* correct the solar constant for the eccentricity of the sun's orbit */
    orbit_correction = 1.0 + 0.033 * cos(doy / 365.0 * 2.0 * M_PI);

    return (1367.0 * orbit_correction * sin_beta);

}


double estimate_clearness(double sw_rad, double So) {
    /* estimate atmospheric transmisivity - the amount of diffuse radiation
    is a function of the amount of haze and/or clouds in the sky. Estimate
    a proxy for this, i.e. the ratio between global solar radiation on a
    horizontal surface at the ground and the extraterrestrial solar
    radiation */
    double tau;

    /* catch possible divide by zero when zenith = 90. */
    if (So <= 0.0)
        tau = 0.0;
    else
        tau = sw_rad / So ;

    if (tau > 1.0) {
        tau = 1.0;
    } else if (tau < 0.0) {
        tau = 0.0;
    }

    return (tau);
}
