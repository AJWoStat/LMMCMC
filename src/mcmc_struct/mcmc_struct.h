#ifndef MCMC_STRUCT_H
#define MCMC_STRUCT_H

#include <stdlib.h>
#include <math.h>
#include "../model_struct/model_struct.h"
#include "../hash_table_struct/hash_table_struct.h"
#include "../mcmc_proposal_functions/mcmc_proposal_functions.h"
// #include "model_struct.h"
// #include "hash_table_struct.h"
// #include "mcmc_proposal_functions.h"

struct mcmc_struct{
  model_t * mod_curr;
  model_t * mod_prop;
  double cumsum_proposal_probabilities[3]; //note the hard-coding of 3 proposals here
  double log_prob_forward;
  double log_prob_backward;
  int num_burn_in;
  int num_thin;
  int num_draws;
  int bitrep_length;
  int next_mcmc_id;
  int draws[];
};

extern struct mcmc_struct * mcmc;

void mcmc_struct_constructor(model_t * model_start, int num_burn_in, int num_thin, int num_draws, double proposal_probabilities[3]);
void mcmc_struct_destructor();

void mcmc_draw(hash_table_t ** ht, int * hash_table_location_curr, int * linked_list_location_curr);
void mcmc_all_draws(hash_table_t ** ht);
int get_mcmc_id_and_increase_count(hash_table_t ** ht, int hash_table_location_value, int linked_list_location_value);



#endif
