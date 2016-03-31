#ifndef READ_MET_H
#define READ_MET_H

#include <stdio.h>
#include <stdlib.h>

#include "gday.h"
#include "utilities.h"

void    read_daily_met_data(char **, control *, met_arrays *);
void    read_subdaily_met_data(char **, control *, met_arrays *);


#endif /* READ_MET_H */
