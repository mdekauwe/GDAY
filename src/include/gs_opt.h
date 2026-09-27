#ifndef GS_OPT_H
#define GS_OPT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "gday.h"
#include "structures.h"
#include "constants.h"
#include "utilities.h"
#include "photosynthesis.h"

/*
** A big leaf (canopy) that supplies its own A(Ci) and E(gs), for the daily
** (MATE) model: gs_opt_canopy() optimises Ci with the same profit and
** hydraulics as the sub-daily leaves.
*/
typedef struct {
    double (*assim)(double ci, void *ctx);   /* A (umol m-2 s-1) at Ci */
    double (*trans)(double gsc, void *ctx);  /* E (mmol m-2 s-1) at gs for
                                                CO2 (mol m-2 s-1) */
    void   *ctx;
    double gamma_star;   /* umol mol-1, lower end of the Ci search */
    double ca;           /* umol mol-1, upper end */
    double psi_rz;       /* root zone water potential (MPa) */
    double kmax;         /* plant conductance, same area basis as E */
    double k_soil;       /* soil to root conductance ahead of the plant,
                            same basis, < 0 none */
    double e_scale;      /* E costed = e_scale x trans(gs) */
    double e_max;        /* trans(gs) can't exceed the soil's supply (mmol
                            m-2 s-1), < 0 no limit */
} gs_opt_canopy_in;

void   gs_opt_leaf(control *, canopy_wk *, met *, params *, state *);
int    gs_opt_canopy(control *, params *, const gs_opt_canopy_in *, double *,
                     double *, double *, double *);
double gs_opt_psi_leaf(control *, canopy_wk *, params *, state *, double,
                       double *);
double gs_opt_beta(control *, canopy_wk *, met *, params *, state *);
void   weibull_params(params *, double *, double *);

#endif /* GS_OPT_H */
