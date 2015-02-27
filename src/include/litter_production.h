#ifndef LITTER_H
#define LITTER_H

#include "gday.h"

/* litter stuff */
void   daily_grazing_calc(double , params *, fluxes *, state *);
void   annual_grazing_calc(params *, fluxes *, state *);
float  decay_in_dry_soils(double, double, params *, state *);
void   calculate_litterfall(control *, fluxes *, params *, state *, int,
                            double *, double *);

#endif /* LITTER */
