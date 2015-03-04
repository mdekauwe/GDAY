#ifndef OPTROOT_H
#define OPTROOT_H

#include "gday.h"

void calc_opt_root_depth(double, double, double, double, double, double,
                         double *, double *, double *);
double estimate_max_root_depth(double, double, double, double);
double rtot_wrapper(double, double, double, double);
double rtot(double, double, double);
double rtot_derivative(double, double, double, double);
double calculate_root_mass_above_depth(double, double, double, double, double);
double calc_plant_nuptake(double, double, double, double);
double calc_umax(double, double, double);
double calc_net_n_uptake(double, double, double, double, double);
double newton(double (*)(double, double, double, double),
              double (*)(double, double, double, double), double, double,
              double, double);
#endif /*  OPTROOT_H */
