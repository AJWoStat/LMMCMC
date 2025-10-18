#ifndef HASH_KEY_COMPUTER_H
#define HASH_KEY_COMPUTER_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "R_ext/Print.h"

struct hash_key_parameters_struct{
  int shift_right_xor[4];
  uint64_t times_replace[3];
  uint64_t vec_add;
  int vec_mult_length;
  uint64_t vec_mult[];
};

static struct hash_key_parameters_struct * hash_key_par;

void hash_key_parameters_constructor(const int vec_mult_length);
void hash_key_parameters_destructor();
void hash_key_parameters_print();

uint32_t triple32inc(uint32_t u);
uint32_t hash_key_computer(const uint32_t * u, const int u_length);

#endif