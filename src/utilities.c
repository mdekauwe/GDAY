
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


double get_diffuse_frac(int doy, double par) {
    /*
        For the moment, I am only going to implement Spitters, so this is a bit
        of a useless wrapper function.

    */

    return spitters(doy, par);
}

double spitters(int doy, double par) {

    /*

    Spitters algorithm to estimate the diffuse component from the measured
    irradiance.

    Parameters:
    ----------
    doy : int
        day of year
    zenith : double
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

    double sw_rad;
    double diffuse = -99.;



    return (diffuse);

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

    double zen_angle, dec, et, lc, t0, h, rlat;

    dec = calculate_solar_declination(doy);

    /* dec = self.calculate_solar_declination(doy)
    et = self.calculate_eqn_of_time(doy)
    lc = self.calculate_longitudal_correction(lon)
    t0 = 12.0 - lc - et # time of solar noon (h)

    # solar hour
    h = 15.0 * (hour - t0) * self.deg2rad

    rlat = lat * self.deg2rad
    cos_zen = (np.sin(rlat) * np.sin(dec) + np.cos(rlat) * np.cos(dec) *
                np.cos(h))
    cos_zen = np.where(cos_zen < 0.0, 0.0, cos_zen)

    zen_angle = np.arccos(cos_zen) * self.rad2deg
    zen_angle = np.ma.masked_array(zen_angle, mask=cos_zen < 0.0001)


    return cos_zen, zen_angle
    */

    return (zen_angle);
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

    sindec = -sin(23.5 * DEG2RAD) * cos(2.0 * PI * (doy + 10.0) / 365.);

    return (ARCSIN(sindec));

}
