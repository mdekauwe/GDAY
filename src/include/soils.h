#ifndef SOILS_H
#define SOILS_H

#include "gday.h"
#include "utilities.h"
#include "constants.h"
#include "water_balance.h"

double calc_soil_temp_factor(double);
void   calculate_csoil_flows(control *, fluxes *, params *, state *,
                             double, int);
void   calculate_decay_rates(fluxes *, params *, state *);
void   flux_from_grazers(control *c, fluxes *, params *);
double calc_ligin_nratio_leaves(control *c, fluxes *f, params *);
double calc_ligin_nratio_fine_roots(control *c, fluxes *f, params *);
double ratio_of_litternc_to_live_rootnc(control *, fluxes *, params *);
double metafract(double);
void   partition_plant_litter(fluxes *, params *);
double ratio_of_litternc_to_live_leafnc(control *, fluxes *, params *);
void   cfluxes_from_structural_pool(fluxes *, params *, state *);
void   cfluxes_from_metabolic_pool(fluxes *f, params *, state *);
void   cfluxes_from_active_pool(fluxes *, params *, state *, double);
void   cfluxes_from_slow_pool(fluxes *, params *, state *);
void   cfluxes_from_passive_pool(fluxes *, params *, state *);
void   calculate_soil_respiration(control *, fluxes *, params *, state *);
void   calculate_cpools(fluxes *, state *);
void   precision_control_soil_c(fluxes *, state *);



/* N stuff */
void   calculate_nsoil_flows(control *, fluxes *, params *, state *, int);
void   grazer_inputs(control *, fluxes *, params *);
void   n_inputs_from_plant_litter(fluxes *, params *, double *, double *);
void   partition_plant_litter_n(control *, fluxes *, params *, double, double);
void   nfluxes_from_structural_pools(fluxes *, params *, state *);
void   nfluxes_from_metabolic_pool(fluxes *, params *, state *);
void   nfluxes_from_active_pool(fluxes *, params *, state *, double);
void   nfluxes_from_slow_pool(fluxes *, params *, state *s);
void   nfluxes_from_passive_pool(fluxes *, params *, state *);
void   calculate_n_mineralisation(fluxes *);
void   calculate_n_immobilisation(fluxes *, params *, state *, double *,
                                  double *, double *, double *);
void   calc_n_net_mineralisation(fluxes *);
double calculate_nc_slope(params *, double, double);
void   calculate_npools(control *c, fluxes *, params *, state *, double,
                        double, double);
double nc_limit(fluxes *, double, double, double, double);
double nc_flux(double, double, double);
void   precision_control_soil_n(fluxes *, state *);

/* P stuff */
void   calculate_psoil_flows(control *, fluxes *, params *, state *, int);
void   grazer_inputs_p(control *, fluxes *, params *);
void   p_inputs_from_plant_litter(fluxes *, params *, double *, double *);
void   partition_plant_litter_p(control *, fluxes *, params *, double, double);
void   pfluxes_from_structural_pools(fluxes *, params *, state *);
void   pfluxes_from_metabolic_pool(fluxes *, params *, state *);
void   pfluxes_from_active_pool(fluxes *, params *, state *, double);
void   pfluxes_from_slow_pool(fluxes *, params *, state *s);
void   pfluxes_from_passive_pool(fluxes *, params *, state *);
void   calculate_p_parent_influx(fluxes *, params *, state *);
void   calculate_p_mineralisation(fluxes *);
void   calculate_p_min_partition(fluxes *, params *, state *);
void   calculate_p_immobilisation(fluxes *, params *, state *, double *,
                                  double *, double *, double *);
void   calculate_p_ssorb_to_sorb(state *, fluxes *, params *, control *);
void   calculate_p_sorb_to_ssorb(state *, fluxes *, params *);
void   calculate_p_ssorb_to_occ(state *, fluxes *, params *); 
void   calc_p_net_mineralisation(fluxes *);
double calculate_pc_slope(params *, double, double);
void   calculate_p_biochemical_mineralisation(fluxes *, params *,state *);
void   calculate_ppools(control *c, fluxes *, params *, state *, double,
                        double, double);
double pc_limit(fluxes *, double, double, double, double);
double pc_flux(double, double, double);
void   precision_control_soil_p(fluxes *, state *);
void   soil_soprtion_parameters(char *, params *);

/* priming/exudation */
void calc_root_exudation_uptake_of_C(fluxes *, params *, state *);
void calc_root_exudation_uptake_of_N(fluxes *, state *);
void calc_root_exudation_uptake_of_P(fluxes *, state *);
void adjust_residence_time_of_slow_pool(fluxes *, params *);

#endif /* SOILS_H */
