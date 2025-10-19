#ifndef MCMC_PROPOSAL_FUNCTIONS_H
#define MCMC_PROPOSAL_FUNCTIONS_H

#include <stdlib.h>
#include <math.h>
#include "Rmath.h"
// #include "gsl/gsl_rng.h"
// #include "gsl/gsl_sf_gamma.h"
#include "../data_store_struct/data_store_struct.h"
#include "../model_struct/model_struct.h"
#include "../hash_table_struct/hash_table_struct.h"
// #include "data_store_struct.h"
// #include "model_struct.h"
// #include "hash_table_struct.h"


// static gsl_rng * rng;
//
// void rng_constructor(unsigned long int seed);
// void rng_destructor();

double draw_unif_continuous(); //just wrapping gsl_unif_pos
int draw_unif_discrete(int upper); //draws from 0 to upper-1, so gives a location in a length upper vector
int draw_discrete(double * cumsum_probs, int length); //draws from 0 to upper-1 with probs given by differences in cumsum probs
int * draw_without_replacement(int size, int upper); //draws size from 0,...,upper-1 without replacement

//modify mod_prop after copy mod_curr into it. return proposal probability.
double simple_random_walk_proposal(const model_t * mod_curr, model_t * mod_prop);
double add_and_remove_proposal(const model_t * mod_curr, model_t * mod_prop);
double random_same_size_bitrep_proposal(const model_t * mod_curr, model_t * mod_prop);

#endif
