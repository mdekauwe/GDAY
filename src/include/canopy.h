#ifndef CANOPY_H
#define CANOPY_H

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "gday.h"
#include "structures.h"
#include "constants.h"
#include "water_balance.h"
#include "water_balance_sub_daily.h"
#include "utilities.h"
#include "radiation.h"
#include "photosynthesis.h"

/* C stuff */
void    initialise_leaf_surface(canopy_wk *, met *);
void    zero_carbon_day_fluxes(fluxes *);
void    zero_hourly_fluxes(canopy_wk *);
void    update_daily_carbon_fluxes(fluxes *, params *, double, double);
void    canopy(canopy_wk *, control *, fluxes *, met_arrays *, met *,
               nrutil *nr, params *, state *);
void    solve_leaf_energy_balance(control *, canopy_wk *, fluxes *, met *,
                                  params *, state *);
void    sum_hourly_carbon_fluxes(canopy_wk *, fluxes *, params *);
void    scale_leaf_to_canopy(control *c, canopy_wk *, state *);
double  calc_leaf_net_rad(params *, state *, double, double, double);
void    calculate_top_of_canopy_leafn(canopy_wk *, params *, state *);
void    calculate_top_of_canopy_leafp(canopy_wk *, params *, state *);
void    calc_leaf_to_canopy_scalar(canopy_wk *, params *);
void    unpack_solar_geometry(canopy_wk *, control *);
/* SPA stuff */
void    check_water_balance(control *, fluxes *, double, double);
double  calc_lwp(fluxes *, state *, double, double);
void    calculate_emax(control *, canopy_wk *, fluxes *, met *, params *, state *);
#endif /* CANOPY_H */
