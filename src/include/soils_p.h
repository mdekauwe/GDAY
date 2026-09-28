#ifndef SOILS_P_H
#define SOILS_P_H

#include "gday.h"
#include "utilities.h"
#include "constants.h"

void   calculate_psoil_flows(control *, fluxes *, params *, state *, int);
double langmuir_labile_p(double, double, double);

#endif /* SOILS_P_H */
