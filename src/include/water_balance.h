#ifndef WATER_BALANCE_H
#define WATER_BALANCE_H

#include "gday.h"
#include "constants.h"
#include "utilities.h"

void   solve_leaf_energy_balance(fluxes *, met *, params *, long, double,
                                 double *, double *, double *, double *);
double  penman_leaf(double, double, double, double, double, double, double,
                    double, double *);
double  calc_radiation_conductance(double);
double  calc_bdn_layer_forced_conduct(double, double, double,double);
double  calc_bdn_layer_free_conduct(double, double, double, double);
double  calc_sat_water_vapour_press(double);
void    calculate_daily_water_balance(control *, fluxes *, met *, params *,
                                      state *, int, double);
void    get_met_data(met *, int , double *, double *, double *, double *,
                    double *, double *, double *, double *, double *, double *,
                    double *, double *, double *, double *, double *);
void    calc_infiltration(fluxes *, params *, state*, double);
double  calc_stomatal_conductance(double, double, double, double, double,
                                 double);
double  calc_radiation(params *, double, double, double);
double  update_water_storage(fluxes *, params *, state *);
void    calc_transpiration_penmon(fluxes *, params *, state *, double, double,
                               double, double, double, double, double);
void    calc_transpiration_penmon_am_pm(params *, state *, double, double, double,
                                       double, double, double, double, double,
                                       double *, double *, double *, double *);
double  calc_soil_evaporation(state*, double, double, double, double);
double  calc_density_of_air(double);
double  calc_latent_heat_of_vapourisation(double);
double  calc_atmos_pressure();
double  calc_pyschrometric_constant(double, double);
double  calc_slope_of_saturation_vapour_pressure_curve(double);
double  canopy_boundary_layer_conductance(params *p, double, double);
void    penman_monteith(double, double, double, double, double, double,
                        double *, double*);
double  calc_sw_modifier(double, double, double);
void    initialise_soil_moisture_parameters(control *, params *);
double *get_soil_fracs(char *);

void    get_soil_params(char *, double *, double *);
void    calc_soil_params(double *, double *, double *,
                        double *, double *, double *);
void    calculate_soil_water_fac(control *, params *, state *);
#endif /* WATER_BALANCE */
