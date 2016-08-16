#ifndef GDAY_H
#define GDAY_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>


#define EPSILON 1E-08
#define DEG2RAD(DEG) (DEG * M_PI / 180.0)
#define RAD2DEG(RAD) (180.0 * RAD / M_PI)
#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#define STRING_LENGTH 2000

#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define CLIP(x) ((x)<0. ? 0. : ((x)>1. ? 1. : (x)))

/* Stomatal conductanct models */
#define MEDLYN 0

/* Photosynthesis models */
#define MATE 0
#define BEWDY 1

/* Respiration models */
#define FIXED 0
#define TEMPERATURE 1
#define BIOMASS 2

/* Allocation models */
#define FIXED 0
#define GRASSES 1
#define ALLOMETRIC 2

/* Ps pathway */
#define C3 0
#define C4 1

/* output time step, where end = the final state */
#define SUBDAILY 0
#define DAILY 1
#define END 2

/* Texture identifiers */
#define SILT 0
#define SAND 1
#define CLAY 2

/* Spinup method */
#define BRUTE 0
#define SAS 1

/* Spinup array index */
#define AF 0
#define AR 1
#define ACR 2
#define AB 3
#define AW 4
#define LF 5
#define LR 6
#define LCR 7
#define LB 8
#define LW 9
#define LSW 10
#define S1 11
#define S2 12
#define KD1 13
#define KD2 14
#define KD3 15
#define KD4 16
#define KD5 17
#define KD6 18
#define KD7 19



#include "structures.h"
#include "initialise_model.h"
#include "simple_moving_average.h"
#include "plant_growth.h"
#include "litter_production.h"
#include "write_output_file.h"
#include "read_param_file.h"
#include "read_met_file.h"
#include "disturbance.h"
#include "phenology.h"
#include "soils.h"
#include "version.h"


void   clparser(int, char **, control *);
void   usage(char **);

void   run_sim(canopy_wk *, control *, fast_spinup *fs, fluxes *,
               met_arrays *, met *, params *p, state *);
void   spin_up_pools(canopy_wk *, control *, fast_spinup *, fluxes *,
                     met_arrays *, met *, params *p, state *);
void   correct_rate_constants(params *, int output);
void   reset_all_n_pools_and_fluxes(fluxes *, state *);
void   zero_stuff(control *, state *);
void   day_end_calculations(control *, params *, state *, int, int);
void   unpack_met_data(control *, fluxes *f, met_arrays *, met *, int, double);
void   zero_fast_spinup_stuff(fast_spinup *);


#endif /* GDAY_H */
