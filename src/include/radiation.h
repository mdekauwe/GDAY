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
void   calculate_solar_geometry(params *, double, double, double *, double *);
double calculate_solar_declination(int, double);
double calculate_eqn_of_time(double);
double get_diffuse_frac(int, double, double);
double spitters(int, double, double);
double calc_extra_terrestrial_irradiance(double, double);
double estimate_clearness(double, double);
void   calculate_absorbed_radiation(params *, state *, double, double, double,
                                    double, double *, double *, double *);
double calculate_solar_noon(double, double);
double calculate_hour_angle(double, double);

#endif /* RADIATION_H */
