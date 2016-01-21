
#include "utilities.h"



int is_leap_year(int yr) {

    if (yr % 400 == 0 || (yr % 100 != 0 && yr % 4 == 0)) {
        return TRUE;
    } else {
        return FALSE;
    }
}



double day_length(int doy, int num_days, double latitude) {

    /*

    Daylength in hours

    Eqns come from Leuning A4, A5 and A6, pg. 1196

    Reference:
    ----------
    Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.

    Parameters:
    -----------
    doy : int
        day of year, 1=jan 1
    yr_days : int
        number of days in a year, 365 or 366
    latitude : float
        latitude [degrees]

    Returns:
    --------
    dayl : float
        daylength [hrs]

    */
    double deg2rad, latr, sindec, a, b;

    deg2rad = M_PI / 180.0;
    latr = latitude * deg2rad;
    sindec = -sin(23.5 * deg2rad) * cos(2.0 * M_PI * (doy + 10.0) / num_days);
    a = sin(latr) * sindec;
    b = cos(latr) * cos(asin(sindec));

    return 12.0 * (1.0 + (2.0 / M_PI) * asin(a / b));
}

void calculate_daylength(int num_days, double latitude, double *dayl) {
    /* wrapper to put the day length into an array */
    int i;
    for (i = 0; i < num_days; i++) {
        dayl[i] = day_length(i+1, num_days, latitude);
    }
    return;
}

void prog_error(const char *reason, const unsigned int line)
{
    fprintf(stderr, "%s, failed at line: %d\n", reason, line);
	exit(EXIT_FAILURE);

    return;
}

bool float_eq(double a, double b) {
    /*
    Are two floats approximately equal...?

    Reference:
    ----------
    D. E. Knuth. The Art of Computer Programming. Sec. 4.2.2 pp. 217-8.
    */
    return fabs(a - b) <= EPSILON * fabs(a);
}



char* rstrip(char* s)
{
    /* Strip whitespace chars off end of given string, in place. Return s. */

    char* p = s + strlen(s);
    while (p > s && isspace(*--p)) *p = '\0';
    return s;
}


char* lskip(char* s)
{
    /* Return pointer to first non-whitespace char in given string. */

    while (*s && isspace(*s)) s++;
    return (char*)s;
}

char* find_char_or_comment(char* s, char c)
{
    /*

    Return pointer to first char c or ';' comment in given string, or
    pointer to null at end of string if neither found. ';' must be
    prefixed by a whitespace character to register as a comment.

    */

    int was_whitespace = 0;
    while (*s && *s != c && !(was_whitespace && *s == ';')) {
        was_whitespace = isspace(*s);
        s++;
    }
    return (char*)s;
}


char *strncpy0(char* dest, char* src, size_t size)
{
    /* Version of strncpy that ensures dest (size bytes) is null-terminated. */

    strncpy(dest, src, size);
    dest[size - 1] = '\0';
    return dest;
}


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
    * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
      biophysics. Pg 168.

    */

    double dec, et, lc, t0, h, gamma, zenith_angle;
    double hour_angle, rlat, rlon, loc_hour_angle, cosz;
    /* need to convert 30 min data, 0-47 to 0-23.5 */

    hod /= 2.0;

    gamma = day_angle(doy);
    dec = calculate_solar_declination(doy, gamma);
    et = calculate_eqn_of_time(doy, gamma);
    lc = (p->longitude - round_to_value(p->longitude, 15.)) / 15.0;
    t0 = 12.0 - lc - et / 60.;
    /*lc = (p->longitude - round_to_value(p->longitude, 15.)) * 4.0;
    t0 = 12.0 - (lc / 60.) - (et / 60.); */
    h = DEG2RAD(15.0 * (hod - t0));
    rlat = DEG2RAD(p->latitude);
    *cos_zen = (sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(h));
    if (*cos_zen < 0.0)
        *cos_zen = 0.0;
    else if (*cos_zen > 1.0)
        *cos_zen = 1.0;

    zenith_angle = RAD2DEG(acos(*cos_zen));
    *elevation = 90.0 - zenith_angle;

    /*printf("%lf %lf %lf %lf\n", hod, zenith_angle, 90.-zenith_angle, cos_zen);*/

    return;
}

double day_angle(int doy) {
    /* Calculation of day angle

    Reference:
    ----------
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.

    Returns:
    ---------
    gamma - day angle in radians.
    */
    double gamma;

    gamma = (2.0 * M_PI * (doy - 1.0)) / 365;

    return (gamma);
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
    * Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    * J. W. Spencer (1971). Fourier series representation of the position of
      the sun.
    */
    double sindec, decl;

    /* declination (radians) */
    decl = 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - \
           0.006758 * cos(2.0 * gamma) + 0.000907 * sin(2.0 * gamma) -\
           0.002697 * cos(3.0 * gamma) + 0.00148 * sin(3.0 * gamma);

    return (decl);

}

double calculate_eqn_of_time(int doy, double gamma) {
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
    et = 0.000075 + 0.001868 * cos(gamma) - 0.032077 * sin(gamma) -\
         0.014615 * cos(2.0 * gamma) - 0.04089 * sin(2.0 * gamma);

    /* radians to minutes */
    et *= 229.18;


    f = 279.575 + 0.9856 * doy;
    A = f * M_PI / 180.0;
    et = (-104.7 * sin(A) + 596.2 * sin(2.0 * A) + 4.3 * sin(3.0 * A) -\
            12.7 * sin(4.0 * A) - 429.3 * cos(A) - 2.0 * cos(2.0 * A) +\
            19.3 * cos(3.0 * A)) / 3600.0;


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

double round_to_value(double number, double roundto) {
    return (round(number / roundto) * roundto);
}
