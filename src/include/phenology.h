#ifndef PHEN_H
#define PHEN_H

#include "gday.h"

void    phenology(control *, fluxes *, met_arrays *, params *, state *);
void    calculate_leafon_off(control *, met_arrays *, params *, state *, double,
                             double, double, double, int, int *, int *,
                             int *, int *, double);
double  calc_gdd(double);
double  gdd_chill_thresh(double, double, double, double);
double  calc_ncd(double);
double  leaf_drop(double, double, double);
void    calc_ini_grass_pheno_stuff(control *, met_arrays *, int, double *, double *,
                                   double *, double *);
void    calculate_growing_season_fluxes(fluxes *f, state *, int);
void    calculate_days_left_in_growing_season(control *, state *, int, int, int);


#endif /* PHEN_H */
