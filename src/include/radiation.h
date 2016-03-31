#ifndef RADIATION_H
#define RADIATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <ctype.h>
#include "gday.h"
#include "constants.h"
#include "utilities.h"

/* utilities */
double day_angle(int);
void   calculate_solar_geometry(canopy_wk *, params *, double, double);
double calculate_solar_declination(int, double);
double calculate_eqn_of_time(double);
void   get_diffuse_frac(canopy_wk *, int, double);
void   spitters(canopy_wk *, int, double);
double calc_extra_terrestrial_rad(double, double);
double estimate_clearness(double, double);
void   calculate_absorbed_radiation(canopy_wk *, params *, state *, double);
double calculate_solar_noon(double, double);
double calculate_hour_angle(double, double);

#endif /* RADIATION_H */
