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
void    update_daily_carbon_fluxes(fluxes *, params *, double, double);
void    canopy(control *, fluxes *, met *, params *, state *);

void    solve_leaf_energy_balance(control *, fluxes *, met *, params *, state *,
                                  double, double, double, double, double,
                                  double *, double *, double *, double *);
double  calc_radiation_conductance(double);
double  calc_bdn_layer_forced_conduct(double, double, double,double);
double  calc_bdn_layer_free_conduct(double, double, double, double);
void    calculate_top_of_canopy_leafn(params *, state *, double, double,
                                      double *);


#endif /* CANOPY_H */
