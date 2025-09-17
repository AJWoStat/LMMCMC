#include "mcmc_struct.h"

struct mcmc_struct * mcmc;

void mcmc_struct_constructor(model_t * model_start, int num_burn_in, int num_thin, int num_draws, double proposal_probabilities[3]){
  //note the hard-coding of 3 proposals here
  mcmc = malloc(sizeof(struct mcmc_struct) + (num_draws)*sizeof(int));
  mcmc->num_burn_in = num_burn_in;
  mcmc->num_thin = num_thin;
  mcmc->num_draws = num_draws;
  mcmc->cumsum_proposal_probabilities[0] = proposal_probabilities[0];
  for(int i=1; i<3; i++) mcmc->cumsum_proposal_probabilities[i] = proposal_probabilities[i] + mcmc->cumsum_proposal_probabilities[i-1];
  for(int i=0; i<3; i++) mcmc->cumsum_proposal_probabilities[i] /= mcmc->cumsum_proposal_probabilities[2];
  mcmc->mod_curr = model_copy_constructor(model_start);
  mcmc->mod_prop = model_copy_constructor(model_start);
  mcmc->bitrep_length = mcmc->mod_curr->bitrep_length;
  mcmc->log_prob_forward = -log(2.00);
  mcmc->log_prob_backward = mcmc->log_prob_forward;
  mcmc->next_mcmc_id = 0;
  for(int i=0; i<mcmc->num_draws; i++) mcmc->draws[i] = -1;
}

void mcmc_struct_destructor(){
  model_destructor(&(mcmc->mod_curr));
  model_destructor(&(mcmc->mod_prop));
  free(mcmc); mcmc = NULL;
}

void mcmc_draw(hash_table_t ** ht, int * hash_table_location_curr, int * linked_list_location_curr){
  //note the hard-coding of 3 proposals here
  int proposal_draw; 
  if(mcmc->mod_curr->size==data.base_model_size || mcmc->mod_curr->size==data.num_var){
    proposal_draw = 0;
  }else{
    proposal_draw = draw_discrete(mcmc->cumsum_proposal_probabilities, 3); 
  }
  switch(proposal_draw){
  case 0:
    mcmc->log_prob_forward = simple_random_walk_proposal(mcmc->mod_curr, mcmc->mod_prop);
    mcmc->log_prob_backward = mcmc->log_prob_forward;
    break;
  case 1:
    mcmc->log_prob_forward = add_and_remove_proposal(mcmc->mod_curr, mcmc->mod_prop);
    mcmc->log_prob_backward = mcmc->log_prob_forward;
    break;
  case 2:
    mcmc->log_prob_forward = random_same_size_bitrep_proposal(mcmc->mod_curr, mcmc->mod_prop);
    mcmc->log_prob_backward = mcmc->log_prob_forward;
    break;
  }
  
  int hash_table_location_prop_value = 0;
  int linked_list_location_prop_value = 0;
  int * hash_table_location_prop = &hash_table_location_prop_value;
  int * linked_list_location_prop = &linked_list_location_prop_value;
  int resize_flag = hash_table_insertion(ht, mcmc->mod_prop, hash_table_location_prop, linked_list_location_prop);
  if(resize_flag==1){
    resize_flag = hash_table_insertion(ht, mcmc->mod_prop, hash_table_location_prop, linked_list_location_prop);
    resize_flag = hash_table_insertion(ht, mcmc->mod_curr, hash_table_location_curr, linked_list_location_curr);
  }
  double log_u = log(draw_unif_continuous());
  double log_mh_alpha = mcmc->mod_prop->log_BF0 + mcmc->mod_prop->log_prior - (mcmc->mod_curr->log_BF0 + mcmc->mod_curr->log_prior) -(mcmc->log_prob_forward - mcmc->log_prob_backward);
  log_mh_alpha = fmin(log_mh_alpha, 0);
  if(log_u<log_mh_alpha){
    model_copy(mcmc->mod_curr, mcmc->mod_prop);
    (*hash_table_location_curr) = hash_table_location_prop_value;
    (*linked_list_location_curr) = linked_list_location_prop_value;
  }
}

int get_mcmc_id_and_increase_count(hash_table_t ** ht, int hash_table_location_value, int linked_list_location_value){
  int loc = hash_table_location_value*((*ht)->linked_list_length)+linked_list_location_value;
  (*ht)->hash_table[loc].mcmc_count++;
  int out = (*ht)->hash_table[loc].mcmc_id;
  if(out == -1){
    out = mcmc->next_mcmc_id;
    mcmc->next_mcmc_id++;
    (*ht)->hash_table[loc].mcmc_id = out;
  }
  return(out);
}

void mcmc_all_draws(hash_table_t ** ht){
  int hash_table_location_value = 0;
  int linked_list_location_value = 0;
  int * hash_table_location = &hash_table_location_value;
  int * linked_list_location = &linked_list_location_value;
  int resize_flag = 1;
  int mcmc_id=-1;
  while(resize_flag==1){
    resize_flag = hash_table_insertion(ht, mcmc->mod_curr, hash_table_location, linked_list_location);
  }
  for(int i=0; i<mcmc->num_burn_in; i++){
    mcmc_draw(ht, hash_table_location, linked_list_location);
  }
  for(int i=0; i<mcmc->num_draws; i++){
    for(int j=0; j<mcmc->num_thin; j++){
      mcmc_draw(ht, hash_table_location, linked_list_location);
    }
    mcmc_id = get_mcmc_id_and_increase_count(ht, hash_table_location_value, linked_list_location_value);
    mcmc->draws[i] = mcmc_id;
  }
}  