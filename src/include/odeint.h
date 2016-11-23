#ifndef ODEINT_H
#define ODEINT_H

#include "gday.h"
//#include "rkqs.h"

void odeint(double [], int, double, double, double,
            double, double, int *, int *,
            double, double, double, double, double, nrutil *,
	        void (*)(double, double [], double [], double, double,
                     double, double, double),
	        void (*)(double [], double [], int, double *, double, double,
					 double [], double *, double *, double, double, double,
                     double, double, nrutil *,
					 void (*)(double, double [], double [], double, double,
 						      double, double, double)));


#endif /* ODEINT_H */
