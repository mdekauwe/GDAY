#ifndef DISTURBANCE_H
#define DISTURBANCE_H

#include "gday.h"


void figure_out_years_with_disturbances(control *, met *, int, int **, int **);
int  time_till_next_disturbance();
int  check_for_fire(int, int, int);
void fire(control *, params *, state *);
void hurricane(fluxes *, params *, state *);

#endif /* DISTURBANCE_H */
