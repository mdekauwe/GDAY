#ifndef ODEINT_H
#define ODEINT_H

#include "gday.h"


void odeint(double [], int, double, double, double, double, double, int *,
            int *, double, double, double,
	        void (*derivs)(double, double [], double [], double, double,
						   double),
	        void (*rkqs)(double [], double [], int, double *, double, double,
						 double [], double *, double *, double, double, double,
						 void (*)(double, double [], double [], double, double,
 						          double)));



/*double zbrent(double (*func) (double, double, double, double, double, double),
              double, double, double, double, double, double, double, double);*/

#endif /* ODEINT_H */
