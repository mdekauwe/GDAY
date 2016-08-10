#ifndef RKQS_H
#define RKQS_H

#include "gday.h"


void (*rkqs)(float [], float [], int, float *, float, float,
             float [], float *, float *, double, double, double,
             void (*)(float, float *, float *,
                      double, double, double));



#endif /* RKQS_H */
