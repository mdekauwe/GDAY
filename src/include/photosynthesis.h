#ifndef PHOTOSYNTHESIS_H
#define PHOTOSYNTHESIS_H

#include "gday.h"
#include "constants.h"
#include "utilities.h"

/* Sub-daily funcs */
void   photosynthesis_C3(control *, canopy_wk *, met *m, params *, state *);
void   photosynthesis_C3_emax(control *, canopy_wk *, met *m, params *,
                              state *);
double calc_co2_compensation_point(params *, double);
double calculate_michaelis_menten(params *, double);
void   calculate_jmaxt_vcmaxt(control *, canopy_wk *, params *, state *,
                              double, double *, double *);
double arrhenius(double, double, double, double);
double peaked_arrhenius(double, double, double, double, double, double);
double calc_leaf_day_respiration(double, double);
double quad(double, double, double, bool, int *);


/* Daily funcs */
void   mate_C3_photosynthesis(control *, fluxes *, met *, params *,
                              state *, double, double, double);

double  calculate_top_of_canopy_n(params *, state *, double);
double  calculate_top_of_canopy_p(params *, state *, double);
double  calculate_co2_compensation_point(params *, double, double);
double  arrh(double, double, double, double);
double  peaked_arrh(double, double, double, double, double, double);
double  calculate_michaelis_menten_parameter(params *, double, double);
void    calculate_jmax_and_vcmax(control *, params *, state *, double, double,
                                 double *, double *, double);
void    calculate_jmax_and_vcmax_with_p(control *, params *, state *, double,
                                        double, double, double *, double *,
                                        double);
void    adj_for_low_temp(double *, double);
double  calculate_ci(control *, params *, state *, double, double);
double  calculate_quantum_efficiency(params *, double ci, double);
double  assim(double, double, double, double);
double  assim_p(double);
double  epsilon(params *, double, double, double, double);


/* C4 additional prototypes */
void   mate_C4_photosynthesis(control *, fluxes *, met *, params *,
                              state *, double, double, double);
void   calculate_vcmax_parameter(params *, state *s, double, double,
                                 double *,
                                 double *, double);
double calc_respiration(double, double);
double quadratic(double, double, double);

#endif /* PHOTOSYNTHESIS */
