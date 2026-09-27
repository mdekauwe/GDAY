#ifndef GS_OPT_H
#define GS_OPT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gday.h"
#include "structures.h"
#include "constants.h"
#include "utilities.h"
#include "photosynthesis.h"

void   gs_opt_leaf(control *, canopy_wk *, met *, params *, state *);
double gs_opt_psi_leaf(control *, canopy_wk *, params *, state *, double,
                       double *);
double gs_opt_beta(control *, canopy_wk *, met *, params *, state *);
void   weibull_params(params *, double *, double *);

#endif /* GS_OPT_H */
