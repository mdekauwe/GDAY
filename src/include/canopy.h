#ifndef CANOPY_H
#define CANOPY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


#include "gday.h"
#include "constants.h"
#include "water_balance.h"
#include "utilities.h"
#include "photosynthesis.h"

/* C stuff */
void    zero_carbon_day_fluxes(fluxes *);
void    update_daily_carbon_fluxes(fluxes *, params *, double);
void    canopy(control *, fluxes *, met *, params *, state *, int);
void    calculate_absorbed_radiation(params *, state *, double, double, double,
                                     double *);

void    solve_leaf_energy_balance(fluxes *, met *, params *, state *, long,
                                  double, double, double, double, double *,
                                  double *, double *, double *);
double  calc_radiation_conductance(double);
double  calc_bdn_layer_forced_conduct(double, double, double,double);
double  calc_bdn_layer_free_conduct(double, double, double, double);



#endif /* CANOPY_H */
