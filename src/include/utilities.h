#ifndef UTILITIES_H
#define UTILITIES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include "gday.h"
#include "constants.h"




/* utilities */
double calculate_zenith_angle(params *, double, double);
double calculate_solar_declination(int);
double calculate_eqn_of_time(int);
double calculate_longitudal_correction(double);
double get_diffuse_frac(int, double, double);
double spitters(int, double, double);
double calc_extra_terrestrial_irradiance(double, double);
double estimate_clearness(double, double);
double day_length(int, int, double);
void   calculate_daylength(int, double, double *);
int    is_leap_year(int);
void   prog_error(const char *, const unsigned int);
bool   float_eq(double, double);

char   *rstrip(char *);
char   *lskip(char *);
char   *find_char_or_comment(char*, char);
char   *strncpy0(char*, char*, size_t);
char   *strip_first_and_last_character(char);

#endif /* UTILITIES_H */
