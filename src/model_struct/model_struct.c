#include "model_struct.h"

void model_alloc_coef_estimates_and_vcov(model_t * mod, const int include_coef_estimates, const int include_coef_vcov){
  // just for use in constructor, bitrep constructor, and intarray constructor
  if(mod->size>0){
    if(include_coef_estimates!=0){
      mod->coef = (double*)calloc(mod->size, sizeof(double));
      mod->coef_sd = (double*)calloc(mod->size, sizeof(double));
      for(int i=0; i<mod->size; i++){mod->coef[i] = 0.0; mod->coef_sd[i] = 0.0;}
      if(include_coef_vcov!=0){
        mod->coef_vcov_upper = (double*)calloc((mod->size<=1 ? 1 : ((mod->size*(mod->size - 1))/2)), sizeof(double));
        for(int i=0; i< (mod->size<=1 ? 1 : ((mod->size*(mod->size - 1))/2)); i++) mod->coef_vcov_upper[i] = 0.0;
      }else mod->coef_vcov_upper = NULL;
    }
  }else{
    if(include_coef_estimates!=0){
      mod->coef = (double*)calloc(1, sizeof(double));
      mod->coef_sd = (double*)calloc(1, sizeof(double));
      mod->coef[0] = 0.0; mod->coef_sd[0] = 0.0;
      if(include_coef_vcov!=0){
        mod->coef_vcov_upper = (double*)calloc(1, sizeof(double));  
        mod->coef_vcov_upper[0] = 0.0;
      }else mod->coef_vcov_upper = NULL;
    }else {mod->coef = NULL; mod->coef_sd = NULL; mod->coef_vcov_upper = NULL;}
  }
  return;
}
void model_realloc_coef_estimates_and_vcov(model_t * mod, const int include_coef_estimates, const int include_coef_vcov){
  // for use in copy and bitrep_replacement and bitflips
  free(mod->coef); mod->coef = NULL;
  free(mod->coef_sd); mod->coef_sd = NULL;
  free(mod->coef_vcov_upper); mod->coef_vcov_upper = NULL;
  model_alloc_coef_estimates_and_vcov(mod, include_coef_estimates, include_coef_vcov);
  return;
}

model_t * model_constructor(const int bitrep_length, const int include_coef_estimates, const int include_coef_vcov){
  model_t * mod = malloc(sizeof(model_t) + bitrep_length*sizeof(uint32_t));
  mod->hash_key = 0;
  mod->size = 0;
  mod->rank = 0;
  mod->Rsq = 0.00;
  mod->log_prior = 0.00;
  mod->log_BF0 = 0.00;
  mod->residual_sd = 0.00;
  mod->shrinkage_factor = log(3.00)/2.00;
  mod->bitrep_length = bitrep_length;
  for(int i=0; i<mod->bitrep_length; i++) mod->bitrep[i]=0;
  model_alloc_coef_estimates_and_vcov(mod, include_coef_estimates, include_coef_vcov);
  return(mod);
}

void model_copy(model_t * mod, const model_t * mod_in){
  if(mod_in == NULL){model_destructor(&mod); return;}
  if(mod->bitrep_length!=mod_in->bitrep_length){
    Rprintf("Bad model_struct copy: bitrep lengths differ. No copy occurred.\n");
    return;
  }
  for(int i=0; i<mod->bitrep_length; i++) mod->bitrep[i] = mod_in->bitrep[i];
  mod->hash_key = mod_in->hash_key;
  mod->size = mod_in->size;
  mod->rank = mod_in->rank;
  mod->Rsq = mod_in->Rsq;
  mod->log_prior = mod_in->log_prior;
  mod->log_BF0 = mod_in->log_BF0;
  mod->residual_sd = mod_in->residual_sd;
  mod->shrinkage_factor = mod_in->shrinkage_factor;
  model_realloc_coef_estimates_and_vcov(mod, mod_in->coef == NULL ? 0 : 1, mod_in->coef_vcov_upper == NULL ? 0 : 1);
  if(mod_in->coef != NULL){
    for(int i=0; i<mod->size; i++){mod->coef[i] = mod_in->coef[i]; mod->coef_sd[i] = mod_in->coef_sd[i];}
    if(mod_in->coef_vcov_upper != NULL) for(int i=0; i< (mod->size<=1 ? 1 : ((mod->size*(mod->size - 1))/2)); i++) mod->coef_vcov_upper[i] = mod_in->coef_vcov_upper[i];
  }
  return;
}
 
model_t * model_copy_constructor(const model_t * mod_in){
  if(mod_in == NULL) return((model_t*)NULL);
  model_t * mod = model_constructor(mod_in->bitrep_length, 0, 0);
  model_copy(mod, mod_in);
  return(mod);
}

model_t * model_bitrep_constructor(const uint32_t * bitrep_in, const int bitrep_length_in, const int include_coef_estimates, const int include_coef_vcov){
  model_t * mod = model_constructor(bitrep_length_in, 0, 0);
  for(int i=0; i<mod->bitrep_length; i++) mod->bitrep[i] = bitrep_in[i];
  mod->size = uint32_vec_count_bits(bitrep_in, bitrep_length_in);
  model_alloc_coef_estimates_and_vcov(mod, include_coef_estimates, include_coef_vcov);
  mod->hash_key = hash_key_computer(mod->bitrep, mod->bitrep_length);
  return(mod);
}

model_t * model_intarray_constructor(const int * intarray_in, const int intarray_length_in, const int bitrep_length_in, const int include_coef_estimates, const int include_coef_vcov){
  model_t * mod = model_constructor(bitrep_length_in, 0, 0);
  mod->size = intarray_length_in;
  model_alloc_coef_estimates_and_vcov(mod, include_coef_estimates, include_coef_vcov);
  uint32_t * bitrep_in = int_to_uint32_t(intarray_in, intarray_length_in, bitrep_length_in);
  for(int i=0; i<mod->bitrep_length; i++) mod->bitrep[i] = bitrep_in[i];
  free(bitrep_in); bitrep_in = NULL;
  mod->hash_key = hash_key_computer(mod->bitrep, mod->bitrep_length);
  return(mod);
}

void model_destructor(model_t ** mod){
  if(*mod != NULL){
    if((*mod)->coef != NULL){free((*mod)->coef); (*mod)->coef = NULL;}
    if((*mod)->coef_sd != NULL){free((*mod)->coef_sd); (*mod)->coef_sd = NULL;}
    if((*mod)->coef_vcov_upper != NULL){free((*mod)->coef_vcov_upper); (*mod)->coef_vcov_upper = NULL;}
    free(*mod); *mod = NULL;
  }
  return;
}

int model_is_equal(const model_t * mod_1, const model_t * mod_2, int only_bitrep){
  if(mod_1==NULL && mod_2!=NULL) return(0);
  if(mod_1!=NULL && mod_2==NULL) return(0);
  if(mod_1==NULL && mod_2==NULL) return(1);
  if(mod_1->bitrep_length != mod_2->bitrep_length) return(0);
  for(int i=0; i<mod_1->bitrep_length; i++) if(mod_1->bitrep[i] != mod_2->bitrep[i]) return(0);
  if(only_bitrep == 0){
    if(mod_1->size != mod_2->size) return(0);
    if(mod_1->rank != mod_2->rank) return(0);
    if(mod_1->hash_key != mod_2->hash_key) return(0);
    if(double_is_equal(mod_1->Rsq, mod_2->Rsq)==0) return(0);
    if(double_is_equal(mod_1->log_prior, mod_2->log_prior)==0) return(0);
    if(double_is_equal(mod_1->log_BF0, mod_2->log_BF0)==0) return(0);
    if(double_is_equal(mod_1->residual_sd, mod_2->residual_sd)==0) return(0);
    if(double_is_equal(mod_1->shrinkage_factor, mod_2->shrinkage_factor)==0) return(0);
    if(mod_1->coef == NULL && mod_2->coef != NULL) return(0);
    if(mod_1->coef != NULL && mod_2->coef == NULL) return(0);
    if(mod_1->coef_sd == NULL && mod_2->coef_sd != NULL) return(0);
    if(mod_1->coef_sd != NULL && mod_2->coef_sd == NULL) return(0);
    if(mod_1->coef_vcov_upper == NULL && mod_2->coef_vcov_upper != NULL) return(0);
    if(mod_1->coef_vcov_upper != NULL && mod_2->coef_vcov_upper == NULL) return(0);
    if(mod_1->coef != NULL) for(int i=0; i<mod_1->size; i++) if(double_is_equal(mod_1->coef[i],mod_2->coef[i])==0) return(0);
    if(mod_1->coef_sd != NULL) for(int i=0; i<mod_1->size; i++) if(double_is_equal(mod_1->coef_sd[i], mod_2->coef_sd[i])==0) return(0);
    if(mod_1->coef_vcov_upper != NULL) for(int i=0; i<(mod_1->size<=1 ? 1: (mod_1->size*(mod_1->size-1))/2); i++) if(double_is_equal(mod_1->coef_vcov_upper[i], mod_2->coef_vcov_upper[i])==0) return(0);
  }
  return(1);
}

void model_bitrep_replacement(model_t * mod, const uint32_t * bitrep, const int bitrep_length){
  //resets model information after bitrep replacement
  if(mod->bitrep_length != bitrep_length){
    Rprintf("Warning: mod->bitrep can't be replaced with a bitrep of a different length (the attempted replacement was ignored).");
    return;
  }
  mod->size = uint32_vec_count_bits(bitrep, bitrep_length);
  for(int i=0; i<bitrep_length; i++) mod->bitrep[i] = bitrep[i];
  model_realloc_coef_estimates_and_vcov(mod, mod->coef==NULL ? 0 : 1, mod->coef_vcov_upper==NULL ? 0 : 1);
  mod->rank = 0;
  mod->Rsq = 0.00;
  mod->log_prior = 0.00;
  mod->log_BF0 = 0.00;
  mod->residual_sd = 0.00;
  mod->shrinkage_factor = log(3.00)/2.00;
  mod->hash_key = hash_key_computer(mod->bitrep, mod->bitrep_length);
  return;
}

void model_bit_flips(model_t * mod, const int * bits, const int bits_length){
  //resets model information after bit flips
  div_t bucket_and_location;
  int bucket, location;
  for(int i=0; i<bits_length; i++){
    bucket_and_location = div(bits[i], 32);
    bucket = bucket_and_location.quot;
    location = bucket_and_location.rem;
    if(bucket > mod->bitrep_length){
      Rprintf("Warning: attempted to flip a bit that was not in mod->bitrep (attempted bit flip was ignored).");
    }else{
      mod->bitrep[bucket] ^= powers_of_2_uint32[location];
    }
  }
  mod->size = uint32_vec_count_bits(mod->bitrep, mod->bitrep_length);
  model_realloc_coef_estimates_and_vcov(mod, mod->coef==NULL ? 0 : 1, mod->coef_vcov_upper==NULL ? 0 : 1);
  mod->rank = 0;
  mod->Rsq = 0.00;
  mod->log_prior = 0.00;
  mod->log_BF0 = 0.00;
  mod->residual_sd = 0.00;
  mod->shrinkage_factor = log(3.00)/2.00;
  mod->hash_key = hash_key_computer(mod->bitrep, mod->bitrep_length);
  return;
}

void model_print(const model_t * mod){
  Rprintf("bitrep_length: %d\n", mod->bitrep_length);
  Rprintf("bitrep:\n");
  for(int i=0; i<mod->bitrep_length; i++ ) Rprintf("%u ", mod->bitrep[i]);
  Rprintf("\n");
  Rprintf("size: %d\n", mod->size);
  Rprintf("hash_key: %u\n", mod->hash_key);
  Rprintf("Rsq: %f\n", mod->Rsq);
  Rprintf("residual_sd: %f\n", mod->residual_sd);
  Rprintf("shrinkage_factor: %f\n", mod->shrinkage_factor);
  Rprintf("log_prior: %f\n", mod->log_prior);
  Rprintf("log_BF0: %f\n", mod->log_BF0);
  if(mod->coef != NULL){
    Rprintf("coef: "); for(int i=0; i<mod->size; i++) Rprintf("%f ", mod->coef[i]); Rprintf("\n");
    Rprintf("coef_sd: "); for(int i=0; i<mod->size; i++) Rprintf("%f ", mod->coef_sd[i]); Rprintf("\n");
    if(mod->coef_vcov_upper != NULL){
      Rprintf("coef_vcov_upper:\n");
      for(int i=0; i<mod->size; i++){
        for(int j=0; j<i; j++) Rprintf("%f ", 0.0);
        for(int j=i; j<mod->size; j++) Rprintf("%f ", mod->coef_vcov_upper[(j*(j+1))/2+i]);
        Rprintf("\n");
      }
      Rprintf("\n");
    }
  }
}
