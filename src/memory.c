


#include "memory.h"

/* Allocate memory for all of the various structures */

void *allocate_struct(size_t struct_size, const unsigned int line)
{
    void *ptr = NULL;

    if (!(ptr = malloc(struct_size))){
        fprintf(stderr, "Error allocating structure on line: %d\n", line);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

char *allocate_memory_char(int size, const unsigned int line)
{
	char *ptr = NULL;
	if (!(ptr = calloc((size), sizeof(char)))) {
		prog_error("Error allocating enough memory for array", line);
	}
	return (ptr);
}


short *allocate_memory_short(int size, const unsigned int line)
{
	short *ptr = NULL;
	if (!(ptr = calloc((size), sizeof(short)))) {
		prog_error("Error allocating enough memory for array", line);
	}
	return (ptr);
}


long *allocate_memory_long(int size, const unsigned int line)
{
	long *ptr = NULL;
	if (!(ptr = calloc((size), sizeof(long)))) {
		prog_error("Error allocating enough memory for array", line);
	}
	return (ptr);
}


int *allocate_memory_int(int size, const unsigned int line)
{
	int *ptr = NULL;
	if (!(ptr = calloc((size), sizeof(int)))) {
		prog_error("Error allocating enough memory for array", line);
	}
	return (ptr);
}


float *allocate_memory_float(int size, const unsigned int line)
{
	float *ptr = NULL;
	if (!(ptr = calloc((size), sizeof(float)))) {
		prog_error("Error allocating enough memory for array", line);
	}
	return (ptr);
}

double *allocate_memory_double(int size, const unsigned int line)
{
	double *ptr = NULL;
	if (!(ptr = (double *)calloc((size), sizeof(double)))) {
		prog_error("Error allocating enough memory for array", line);
	}
	return (ptr);
}
