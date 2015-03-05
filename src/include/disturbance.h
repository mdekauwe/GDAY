#ifndef DISTURBANCE_H
#define DISTURBANCE_H

#include "gday.h"
#include "constants.h"
#include "utilities.h"


void figure_out_years_with_disturbances(control *, met *, params *, int **,
                                        int *);
int  time_till_next_disturbance();
int  check_for_fire(control *, fluxes *f, params *, state *, int, int *, int);
void fire(control *, fluxes *f, params *, state *);
void hurricane(fluxes *, params *, state *);

#endif /* DISTURBANCE_H */
