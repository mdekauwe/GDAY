#ifndef WATER_BALANCE_H
#define WATER_BALANCE_H

#include "gday.h"
#include "constants.h"
#include "utilities.h"

void    update_water_storage(control *, fluxes *, params *, state *, double,
                             double, double, double *, double *, double *,
                             double *);
double  calc_canopy_evaporation(met *, params *, state *, double);
void    calculate_water_balance(control *, fluxes *, met *, params *,
                              state *, int, double, double, double);
void    zero_water_day_fluxes(fluxes *);
void    update_water_storage_recalwb(control *, fluxes *, params *, state *,
                                     met *);
double  calc_soil_evaporation(met *, params*, state *, double);
void    calc_interception(control *c, met *m, params *, fluxes *, state *,
                          double *, double *, double *);
void    penman_canopy_wrapper(params *, state *, double, double, double, double,
                              double, double, double, double *, double *,
                              double *, double *, double *);
void    penman_leaf_wrapper(met *, params *, state *, double, double,
                            double, double *, double *, double *, double *,
                            double *, double *);
void    penman_monteith(double, double, double, double, double, double, double *,
                        double *, double *, double *);
double  calc_sat_water_vapour_press(double);
void    calculate_daily_water_balance(control *, fluxes *, met *, params *,
                                      state *, int, double);
double  calc_stomatal_conductance(params *, state *, double, double, double);
double  calc_net_radiation(params *, double, double);
double  calc_latent_heat_of_vapourisation(double);
double  calc_pyschrometric_constant(double, double);
double  calc_slope_of_sat_vapour_pressure_curve(double);
void    calc_soil_water_potential(control *, params *, state *);
double  calc_sw_modifier(double, double, double);
void    initialise_soil_moisture_parameters(control *, params *);
double *get_soil_fracs(char *);
double  calc_beta(double, double, double, double, double);
void    get_soil_params(char *, double *, double *);
void    calc_soil_params(double *, double *, double *,
                        double *, double *, double *);
void    calculate_soil_water_fac(control *, params *, state *);
void    sum_hourly_water_fluxes(fluxes *, double, double, double, double,
                                double, double, double, double);
void    update_daily_water_struct(fluxes *, double, double, double, double,
                                  double, double, double);
double  calc_radiation_conductance(double);
double  calc_bdn_layer_forced_conduct(double, double, double,double);
double  calc_bdn_layer_free_conduct(double, double, double, double);

double  canopy_boundary_layer_conduct(params *, double, double, double, double);

/* hydraulics stuff */
void    calc_saxton_parameters(params *, double *);

#endif /* WATER_BALANCE */
