#ifndef READ_MET_H
#define READ_MET_H

#include <stdio.h>
#include <stdlib.h>

#include "gday.h"
#include "utilities.h"

void    read_daily_met_data(char **, control *, met_arrays *);
void    read_subdaily_met_data(char **, control *, met_arrays *);
void    read_daily_met_data_binary(char **, control *, met_arrays *);
void    read_subdaily_met_data_binary(char **, control *, met_arrays *,
                                      params *, state *);
int    rand_int(unsigned int, unsigned int);

void   estimate_dirunal_par(float, float, int, float, float *);
void   estimate_diurnal_vph(float, float, float, float, float *);
void   disaggregate_rainfall(float, float *rain);
void   estimate_diurnal_temp(float, float, float, float *);
float  calc_vpd(float, float);


float  spittersx(int, float, float *);
float  day_anglex(int);
float  calculate_solar_declinationx(int, float);
float  calculate_eqn_of_timex(float);
float  calculate_solar_noonx(float, float);
float  calculate_hour_anglex(float, float);
float  calc_extra_terrestrial_radx(int, float);
float  round_to_valuex(float, float);
void   calculate_solar_geometryx(int, float, float, float *);

#endif /* READ_MET_H */
