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
void    calculate_water_balance_sub_daily(control *, canopy_wk *, fluxes *, met *,
                                          nrutil *, params *, state *, int,
                                          double, double, double, double,
                                          double);
void    setup_hydraulics_arrays(fluxes *, params *, state *);

void    sum_hourly_water_fluxes(fluxes *, double, double, double, double,
                                double, double, double, double, double);

void    calc_saxton_stuff(params *, double *);
double  saxton_field_capacity(double, double, double, double, double, double);
double  calc_soil_conductivity(double, double, double, double);
void    calc_soil_water_potential(fluxes *, params *, state *);
void    calc_soil_root_resistance(control *, fluxes *, params *, state *);
void    calc_water_uptake_per_layer(control *, fluxes *, params *, state *);
void    calc_wetting_layers(fluxes *, params *, state *, double, double);
double  calc_infiltration(fluxes *, params *, state *, double);
void    calc_soil_balance(fluxes *, nrutil *, params *, state *, int );
void    calc_soil_balance_cascading(fluxes *, params *, state *, int);
void    soil_water_store(double, double [], double [], double, double, double,
                         double, double);

void   zero_water_movement(fluxes *, params *);
void   extract_water_from_layers(fluxes *, state *, double, double);
double root_zone_supply(fluxes *, state *);
void   update_soil_water_storage(fluxes *, params *, state *, double *, double *);
double calc_qe_flux(fluxes *, params *, state *, double, double, double, double,
                    double);
double calc_soil_boundary_layer_conductance(double, double);
double  soil_psi_raw(params *, int, double);
double  soil_psi(control *, params *, state *, int, double);
double  soil_theta_at_psi(params *, int, double);
double  soil_conductivity(params *, int, double);
void    setup_soil_hydraulics(control *, params *, double *);
void    setup_soil_retention(control *, params *);
void    calc_soil_balance_richards(fluxes *, params *, state *);

#endif /* WATER_BALANCE_SUBDAILY_H */
