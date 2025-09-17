#ifndef INT_BIT_CONVERTERS_H
#define INT_BIT_CONVERTERS_H

#include <stdint.h>
#include <stdlib.h>

extern const uint32_t powers_of_2_uint32[32];

uint32_t * int_to_uint32_t(const int * x, const int x_length, const int out_length);

int * uint32_t_to_int(const uint32_t * u, const int u_length, int * out_length);

// int * uint32_t_to_int_with_tmp(const uint32_t * u, const int u_length, int * out_length, int ** tmp);

int uint32_vec_count_bits(const uint32_t *u, const int u_length);

#endif
