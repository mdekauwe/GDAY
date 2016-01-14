#ifndef READ_MET_H
#define READ_MET_H

#include <stdio.h>
#include <stdlib.h>

#include "gday.h"
#include "utilities.h"

void    read_daily_met_data(char **, control *, met *);
void    read_subdaily_met_data(char **, control *, met *);
double  get_diffuse_frac(int, double);
double  spitters(int, double);

#endif /* READ_MET_H */
