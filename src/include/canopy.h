#ifndef CANOPY_H
#define CANOPY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


#include "gday.h"
#include "constants.h"
#include "water_balance.h"
#include "utilities.h"
#include "radiation.h"
#include "photosynthesis.h"

/* C stuff */
void    zero_carbon_day_fluxes(fluxes *);
void    zero_hourly_fluxes(double *, double *, double *, double *, double *,
                           double *);
void    update_daily_carbon_fluxes(fluxes *, params *, double, double);
void    canopy(control *, fluxes *, met *, params *, state *, double *,
               double *);

void    solve_leaf_energy_balance(control *, fluxes *, met *, params *,
                                  state *, double, double, double, double,
                                  double *, double *, double *, double *,
                                  double *, double *);
void    sum_hourly_carbon_fluxes(fluxes *, params *, double *, double *,
                                 double *);
double  calc_leaf_net_rad(params *, state *, double, double, double);
double  calculate_top_of_canopy_leafn(params *, state *);
#endif /* CANOPY_H */
