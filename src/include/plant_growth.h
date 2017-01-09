#ifndef PLANT_GROWTH_H
#define PLANT_GROWTH_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>


#include "gday.h"
#include "constants.h"
#include "water_balance.h"
#include "utilities.h"
#include "photosynthesis.h"
#include "optimal_root_model.h"
#include "canopy.h"

/* C stuff */
void    calc_day_growth(canopy_wk *, control *, fluxes *, met_arrays *ma, met *,
                        nrutil *, params *, state *, double, int, double,
                        double);
void    carbon_allocation(control *, fluxes *, params *, state *, double, int);
void    calc_carbon_allocation_fracs(control *c, fluxes *, params *, state *,
                                     double);
double  alloc_goal_seek(double, double, double, double);
void    update_plant_state(control *, fluxes *, params *, state *,
                           double, double, int);
void    precision_control(fluxes *, state *);
void    calculate_cnp_store(control *, fluxes *, state *);
void    calculate_average_alloc_fractions(fluxes *, state *, int );
void    allocate_stored_cnp(fluxes *f, params *p, state *s);
void    carbon_daily_production(control *, fluxes *, met *m, params *, state *,
                                double);
void    calculate_subdaily_production(control *, fluxes *, met *m, params *,
                                      state *, int, double);
void    calculate_cnp_wood_ratios(control *c, params *, state *, double, double,
                                  double, double *, double *, double *,
                                  double *, double *,double *, double *,
                                  double *);

/* N stuff */
int    np_allocation(control *c, fluxes *, params *, state *, double,
                     double, double, double, double,
                     double, double, double, double, double, int);
double calculate_nuptake(control *, params *, state *);

double nitrogen_retrans(control *, fluxes *, params *, state *,
                        double, double, int);


/* P stuff */
double calculate_growth_stress_limitation(params *, state *, control *);
double calculate_puptake(control *, params *, state *, fluxes *);
double phosphorus_retrans(control *, fluxes *, params *, state *,
                          double, double, int);

int cut_back_production(control *, fluxes *, params *, state *, double,
                        double, double, double, double, int);

/* Priming/Exudation stuff */
void   calc_root_exudation(control *c, fluxes *, params *p, state *);

/* hydraulics */
void   initialise_roots(fluxes *, params *, state *);
void   update_roots(control *, params *, state *);
double calc_root_dist(double, double, double, double, double, double);

#endif /* PLANT_GROWTH */
