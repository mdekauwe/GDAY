#ifndef PLANT_GROWTH_H
#define PLANT_GROWTH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


#include "gday.h"
#include "constants.h"
#include "water_balance.h"
#include "utilities.h"



void   calc_day_growth(control *, fluxes *, met *, params *, state *, int,
                        double, int, double, double);


void   calculate_water_balance(control *, fluxes *, met *, params *,
                               state *, int, double);
void   mate_C3_photosynthesis(control *, fluxes *, met *, params *,
                              state *, int, double, double);
void   mate_C4_photosynthesis(control *, fluxes *, met *, params *,
                              state *, int, double, double);

int    nitrogen_allocation(control *c, fluxes *, params *, state *, double,
                           double, double, double, double, double, int);
double calculate_growth_stress_limitation(params *, state *);
void   calc_carbon_allocation_fracs(control *c, fluxes *, params *, state *,
                                    double);
void   carbon_allocation(control *, fluxes *, params *, state *,
                       double, int);
double alloc_goal_seek(double, double, double, double);
void   update_plant_state(control *, fluxes *, params *, state *,
                          double, double, int);
void   precision_control(fluxes *, state *);
void   calculate_cn_store(fluxes *, state *);
void   calculate_average_alloc_fractions(fluxes *, state *, int );
void   allocate_stored_c_and_n(fluxes *f, params *p, state *s);

double calculate_nuptake(control *, params *, state *);
void   carbon_production(control *, fluxes *, met *m, params *, state *, int,
                         double);
double nitrogen_retrans(control *, fluxes *, params *, state *,
                        double, double, int);
void   calculate_ncwood_ratios(control *c, params *, state *, double, double *,
                              double *, double *, double *);
void calc_opt_root_depth(double, double, double, double *, double *, double *);
#endif /* PLANT_GROWTH */
