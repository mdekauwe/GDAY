
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


double get_diffuse_frac(int doy, double zenith_angle, double par) {
    /*
        For the moment, I am only going to implement Spitters, so this is a bit
        of a useless wrapper function.

    */

    return spitters(doy, zenith_angle, par);
}

double spitters(int doy, double zenith_angle, double par) {

    /*

    Spitters algorithm to estimate the diffuse component from the measured
    irradiance.

    Parameters:
    ----------
    doy : int
        day of year
    zenith_angle : double
        sun zenith angle [degrees]
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

    double sw_rad, cos_zenith, S_0, tau, R, K, diffuse_frac;

    sw_rad = par * PAR_2_SW;
    cos_zenith = cos(zenith_angle);
    S_0 = calc_extra_terrestrial_irradiance(doy);
    tau = estimate_clearness(sw_rad, S_0, cos_zenith);

    /* the ratio between diffuse and total Solar irradiance (R), eqn 20 */
    R = 0.847 - 1.61 * cos_zenith + 1.04 * (cos_zenith * cos_zenith);
    K = (1.47 - R) / 1.66;

    /*
        Relation btw diffuse frac and atmospheric transmission for hourly
        radiation values, eqn 20a-d
    */
    if (tau <= 0.22) {
        diffuse_frac = 1.0;
    } else if (tau > 0.222 && tau <= 0.35) {
        diffuse_frac = 1.0 - 6.4 * (tau - 0.22) * (tau - 0.22);
    } else if (tau > 0.35 && tau <= K) {
        diffuse_frac = 1.47 - 1.66 * tau;
    } else if (tau > K) {
        diffuse_frac = R;
    }

    if (diffuse_frac <= 0.0) {
        diffuse_frac = 0.0;
    } else if (diffuse_frac >= 1.0) {
        diffuse_frac = 1.0;
    }


    return (diffuse_frac);

}


double calculate_zenith_angle(params *p, double doy, double hod) {

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

    double zenith_angle, dec, et, lc, t0, h, rlat, cos_zen;

    /* need to convert 30 min data, 0-47 to 0-23.5 */
    hod /= 2.0;

    dec = calculate_solar_declination(doy);
    et = calculate_eqn_of_time(doy);
    lc = calculate_longitudal_correction(p->longitude);
    t0 = 12.0 - lc - et; /* time of solar noon (h) */
    h = DEG2RAD(15.0 * (hod - t0));
    rlat = DEG2RAD(p->latitude);
    cos_zen = (sin(rlat) * sin(dec) + cos(rlat) * cos(dec) * cos(h));
    if (cos_zen < 0.0)
        cos_zen = 0.0;

    zenith_angle = RAD2DEG(acos(cos_zen));

    return (zenith_angle);
}


double calculate_solar_declination(int doy) {
    /*
    Solar Declination Angle

    Arguments:
    ----------
    doy : int
        day of year, 1=jan 1

    Returns:
    --------
    dec: float
        Solar Declination Angle [radians]

    Reference:
    ----------
    Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    */
    double sindec;

    sindec = -sin(DEG2RAD(23.5)) * cos(2.0 * M_PI * (doy + 10.0) / 365.);

    return (asin(sindec));

}

double calculate_eqn_of_time(int doy) {
    /* Equation of time - correction for the difference btw solar time
    and the clock time.

    Arguments:
    ----------
    doy : int
        day of year

    References:
    -----------
    * Campbell, G. S. and Norman, J. M. (1998) Introduction to environmental
      biophysics. Pg 169.
    */
    double f, arg1, arg2, arg3, et;

    f = DEG2RAD(279.575 + 0.9856 * doy);
    arg1 = -104.7 * sin(f) + 596.2 * sin(2.0 * f) + 4.3;
    arg2 = sin(3.0 * f) - 12.7 * sin(4.0 * f) - 429.3;
    arg3 = cos(f) - 2.0 * cos(2.0 * f) + 19.3 * cos(3.0 * f);
    et = arg1 * arg2 * arg3;

    return (et / 3600.0);
}


double calculate_longitudal_correction(double lon) {
    /*
    Calculate the longitudal correction for travelling east and west.
    +4 minutes (+1/15 hour) for every degree east of the standard meridian
    and -4 mins for each degree west.

    Arguments:
    ----------
    lon : double
        day of year

    */
    double merid, Ih, A, SM, lc;

    merid = floor(lon / 15.0) * 15.0;
    if (merid < 0.0)
        merid += 15.0;

    Ih = floor(lon / 15.0);
    A = lon / 15.;
    if (A > (Ih + 0.5)) {
        SM = (Ih + 1.0) * 15.0;
    } else {
        SM = Ih * 15.0;
    }

    /* longitudinal correction */
    lc = (lon - SM) * -4.0 / 60.0;

    return (lc);
}

double calc_extra_terrestrial_irradiance(double doy) {
    /* Solar radiation incident outside the earth's atmosphere, e.g.
    extra-terrestrial radiation. The value varies a little with the earths
    orbit.

    Reference:
    ----------
    Leuning et al (1995) Plant, Cell and Environment, 18, 1183-1200.
    */
    double Sc = 1367.0; /* solar constant */
    return (Sc * (1.0 + 0.033 * cos(2.0 * M_PI * (doy - 10.0) / 365.0)));
}

double estimate_clearness(double sw_rad, double S_0, double cos_zenith) {
    /* estimate atmospheric transmisivity - the amount of diffuse radiation
    is a function of the amount of haze and/or clouds in the sky. Estimate
    a proxy for this, i.e. the ratio between global solar radiation on a
    horizontal surface at the ground and the extraterrestrial solar
    radiation */
    double tau;

    tau = sw_rad / (S_0 * cos_zenith);
    if (tau > 1.0) {
        tau = 1.0;
    } else if (tau < 0.0) {
        tau = 0.0;
    }

    return (tau);
}
