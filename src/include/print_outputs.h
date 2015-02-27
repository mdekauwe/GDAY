#ifndef PRINT_OUTPUTS_H
#define PRINT_OUTPUTS_H


#include "gday.h"
#include "utilities.h"

void open_output_file(control *, char *);
void write_output_header(control *);
void write_daily_outputs(control *, fluxes *, state *, int, int);
int  write_final_state(control *, params *p, state *);
int  ohandler(char *, char *, char *, control *, params *p, state *, int *);
#endif /* PRINT_OUTPUTS_H */
