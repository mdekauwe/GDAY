#ifndef NUMERICAL_LIBS_H
#define NUMERICAL_LIBS_H

#include "gday.h"

double zbrent(double (*func) (double, double, double, double, double, double),
                      double, double, double, double, double, double, double,
                      double);


#endif /* NUMERICAL_LIBS_H */
