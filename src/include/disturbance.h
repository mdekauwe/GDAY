#ifndef DISTURBANCE_H
#define DISTURBANCE_H

#include "gday.h"


int  time_till_next_disturbance();
void hurricane(fluxes *, params *, state *);
void fire(control *, params *, state *);
#endif /* DISTURBANCE_H */
