#ifndef WATER_BALANCE_SUBDAILY_H
#define WATER_BALANCE_SUBDAILY_H

#include "gday.h"
#include "constants.h"
#include "utilities.h"
#include "water_balance.h"
#include "zbrent.h"
#include "nrutil.h"
#include "odeint.h"
//#include "rkck.h"
#include "rkqs.h"



void    initialise_soils_sub_daily(control *, fluxes *, params *, state *);
void    calculate_water_balance_sub_daily(control *, fluxes *, met *,
                                           params *, state *, int, double,
                                           double, double);
void    setup_hydraulics_arrays(fluxes *, params *, state *);


void    sum_hourly_water_fluxes(fluxes *, double, double, double, double,
                                double, double, double, double, double);

void    calc_saxton_stuff(params *, double *);
double  saxton_field_capacity(double, double, double, double, double, double);
double  calc_soil_conductivity(double, double, double, double);
void    calc_soil_water_potential(fluxes *, params *, state *);
void    calc_soil_root_resistance(fluxes *, params *, state *);
void    calc_water_uptake_per_layer(fluxes *, params *, state *);
void    calc_wetting_layers(fluxes *, params *, state *, double, double);
double  calc_infiltration(fluxes *, params *, state *, double);
void    calc_soil_balance(fluxes *, params *, state *, int );
void    soil_water_store(double, double [], double [], double, double, double,
                         double, double);



#endif /* WATER_BALANCE_SUBDAILY_H */
