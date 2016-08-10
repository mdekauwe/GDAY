#ifndef ODEINT_H
#define ODEINT_H

#include "gday.h"
#include "nrutil.h"
#include <math.h>

void odeint(double ystart[], int nvar, double x1, double x2, double eps,
            double h1, double hmin, int *nok, int *nbad,
            void (*derivs)(double, int, double[], double[]),
            void (*rkqs)(double[], double[], int, double*, double,
            double, double[], double*, double*,
            void (*)(double, int, double[], double[])));

#endif /* ODEINT_H */
