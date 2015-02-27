#ifndef MEMORY_H
#define MEMORY_H

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "utilities.h"

/* memory stuff */
void   *allocate_struct(size_t, const unsigned int);
char   *allocate_memory_char(int, const unsigned int);
short  *allocate_memory_short(int, const unsigned int);
long   *allocate_memory_long(int, const unsigned int);
int    *allocate_memory_int(int, const unsigned int);
float  *allocate_memory_float(int, const unsigned int);
double *allocate_memory_double(int, const unsigned int);

#endif /* MEMORY_H */
