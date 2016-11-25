#ifndef INI_MODEL_H
#define INI_MODEL_H

#include "gday.h"

void    initialise_control(control *);
void    initialise_params(params *);
void    initialise_fluxes(fluxes *);
void    initialise_state(state *f);
void    initialise_nrutil(nrutil *);

#endif /* INI_MODEL_H */
