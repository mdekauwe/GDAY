#ifndef RKCK_H
#define RKCK_H

#include "gday.h"


void (*rkck)(float [], float [], int, float, float, float [], float [],
             double, double, double,
	          void (*)(float, float *, float *, double, double, double));


#endif /* RKCK_H */
