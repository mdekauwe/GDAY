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
void   calculate_solar_geometry(canopy_wk *, params *, double, double);
void   get_diffuse_frac(canopy_wk *, int, double);
void   spitters(canopy_wk *, int, double);
void   calculate_absorbed_radiation(canopy_wk *, params *, state *, double, double);
double psi_func(double, double);

#endif /* RADIATION_H */
