#include "mcmc_proposal_functions.h"

// static gsl_rng * rng;
//
// void rng_constructor(unsigned long int seed){
//   rng = gsl_rng_alloc(gsl_rng_mt19937);
//   gsl_rng_set(rng, seed);
//   return;
// }
// void rng_destructor(){
//   gsl_rng_free(rng); rng = NULL;
//   return;
// }

double draw_unif_continuous(){
  return(runif(0, 1));
}

int draw_unif_discrete(int upper){
  double u = runif(0,1);
  return((int)trunc(u*(double)upper));
}

int draw_discrete(double * cumsum_probs, int length){
  double u = runif(0,1);
  for(int out=0; out<length; out++) if(u<cumsum_probs[out]) return(out);
  return(length);
}

int * draw_without_replacement(int size, int upper){
  //does not check for k<p, so be careful in calling this
  int * draws_01 = (int*)calloc(upper,sizeof(int));
  int * draws_int = malloc(size*sizeof(int));
  int raw_draw;
  int count;
  for(int i=0; i<size; i++){
    raw_draw = draw_unif_discrete(upper-i);
    int count = -1;
    for(int j=0; j<upper; j++){
      if(draws_01[j]==0){
        count++;
      }
      if(raw_draw==count){
        draws_01[j]++;
        draws_int[i] = j;
        break;
      }
    }
  }
  free(draws_01); draws_01=NULL;
  return(draws_int);
}

double simple_random_walk_proposal(const model_t * mod_curr, model_t * mod_prop){
  model_copy(mod_prop, mod_curr);
  int bit = data.test_var_indices[draw_unif_discrete(data.num_test_var)];
  model_bit_flips(mod_prop, &bit, 1);
  return(-log((double)data.num_test_var));
  //in theory, this function should have the same forward and backward proposal probabilities
  //so this is just here to account for other possible proposal functions in the future
}

double add_and_remove_proposal(const model_t * mod_curr, model_t * mod_prop){
  if(mod_curr->size == data.num_var || mod_curr->size == data.base_model_size){
    //just to catch if any issues happen so that a bad call does not give an error
    //this function should not be called if either of the above conditions are met
    return(simple_random_walk_proposal(mod_curr, mod_prop));
  }
  model_copy(mod_prop, mod_curr);
  int raw_bit_rm = draw_unif_discrete(mod_curr->size - data.base_model_size);
  int raw_bit_add = draw_unif_discrete(data.num_var - mod_curr->size);
  int * bits = malloc(2*sizeof(int));
  int count_in=-1;
  int count_out=-1;
  int var;
  uint32_t flag;
  div_t bucket_and_location;
  int bucket;
  int location;
  for(int i=0; i<data.num_test_var; i++){
    var = data.test_var_indices[i];
    bucket_and_location = div(var, 32);
    bucket = bucket_and_location.quot;
    location = bucket_and_location.rem;
    flag = mod_curr->bitrep[bucket];
    if(location>0) flag >>= location;
    flag &= (uint32_t) 1;
    if(flag== (uint32_t) 0){
      if(count_out<raw_bit_add){
        count_out++;
        if(count_out==raw_bit_add){
          bits[0] = var;
        }
      }
    }else{
      if(count_in<raw_bit_rm){
        count_in++;
        if(count_in==raw_bit_rm){
          bits[1] = var;
        }
      }
    }
    if(count_in == raw_bit_rm && count_out == raw_bit_add) break;
  }
  model_bit_flips(mod_prop, bits, 2);
  free(bits); bits=NULL;
  return(-log((double)(mod_curr->size - data.base_model_size)) - log((double)(data.num_var-mod_curr->size)));
  //in theory, this function should have the same forward and backward proposal probabilities
  //so this is just here to account for other possible proposal functions in the future
}

double random_same_size_bitrep_proposal(const model_t * mod_curr, model_t * mod_prop){
  if(mod_curr->size == data.num_var || mod_curr->size == data.base_model_size){
    //just to catch if any issues happen so that a bad call does not give an error
    //this function should not be called if either of the above conditions are met
    return(simple_random_walk_proposal(mod_curr, mod_prop));
  }
  model_copy(mod_prop, mod_curr);
  for(int i=0; i<mod_prop->bitrep_length; i++) (mod_prop)->bitrep[i] = data.base_model_bitrep[i];
  int num_to_draw = mod_curr->size - data.base_model_size;
  int * bits_raw = draw_without_replacement(num_to_draw, data.num_test_var);
  int * bits = malloc(num_to_draw*sizeof(int));
  for(int i=0; i<num_to_draw; i++) bits[i] = data.test_var_indices[bits_raw[i]];
  model_bit_flips(mod_prop, bits, num_to_draw);
  free(bits_raw); bits_raw=NULL;
  free(bits); bits=NULL;
  return(-lchoose((double)data.num_test_var, (double)num_to_draw));
  //in theory, this function should have the same forward and backward proposal probabilities
  //so this is just here to account for other possible proposal functions in the future
}

