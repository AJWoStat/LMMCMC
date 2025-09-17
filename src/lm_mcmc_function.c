#include "lm_mcmc_function.h"

// [[register]]
SEXP lm_mcmc_function(SEXP X_in, SEXP y_in, SEXP weights_in, SEXP base_model_indices, //data info
              // SEXP bf_type, SEXP bf_params, SEXP scale_by_rank,  //coefficient prior info - for later
              SEXP msp_family, SEXP msp_params, SEXP trunc, //model space prior info, value of trunc is ignored unless a truncated family is specified
              SEXP include_coef, SEXP include_vcov, //return info (must be 0 for now)
              SEXP start_model_indices, //start model
              SEXP n_burnin, SEXP n_thin, SEXP n_draws, SEXP proposal_probs, //mcmc info
              SEXP hash_table_length, SEXP linked_list_length, SEXP max_load_factor //hash table info
){
  //read in and protect
  SEXP X_in_R = PROTECT(duplicate(X_in));
  SEXP y_in_R = PROTECT(duplicate(y_in));
  SEXP weights_in_R = PROTECT(duplicate(weights_in));
  SEXP base_model_indices_R = PROTECT(duplicate(base_model_indices));
  SEXP msp_family_R = PROTECT(duplicate(msp_family));
  SEXP msp_params_R = PROTECT(duplicate(msp_params));
  SEXP trunc_R = PROTECT(duplicate(trunc));
  SEXP include_coef_R = PROTECT(duplicate(include_coef));
  SEXP include_vcov_R = PROTECT(duplicate(include_vcov));
  SEXP start_model_indices_R = PROTECT(duplicate(start_model_indices));
  SEXP n_burnin_R = PROTECT(duplicate(n_burnin));
  SEXP n_thin_R = PROTECT(duplicate(n_thin));
  SEXP n_draws_R = PROTECT(duplicate(n_draws));
  SEXP proposal_probs_R = PROTECT(duplicate(proposal_probs));
  // SEXP seed_R = PROTECT(duplicate(seed));
  SEXP hash_table_length_R = PROTECT(duplicate(hash_table_length));
  SEXP linked_list_length_R = PROTECT(duplicate(linked_list_length));
  SEXP max_load_factor_R = PROTECT(duplicate(max_load_factor));

  GetRNGstate();

  // getting data store
  int num_obs_in = LENGTH(y_in);
  int num_var_in = LENGTH(X_in)/num_obs_in + 1; //adding one for the intercept
  double * X_in_C = malloc(num_obs_in*num_var_in*sizeof(double));
  for(int i=0; i<num_obs_in; i++) X_in_C[i] = 1.00; //adding the intercept
  for(int j=0, i=num_obs_in; i<(num_obs_in*num_var_in); j++, i++) X_in_C[i] = REAL(X_in_R)[j];//getting the rest
  double * y_in_C = malloc(num_obs_in*sizeof(double));
  for(int i=0; i<num_obs_in; i++) y_in_C[i] = REAL(y_in_R)[i];
  double * weights_in_C = malloc(num_obs_in*sizeof(double));
  for(int i=0; i<num_obs_in; i++) weights_in_C[i] = REAL(weights_in_R)[i];
  int base_model_size = LENGTH(base_model_indices) + 1;
  int * base_model_indices_C = malloc(base_model_size*sizeof(int));
  base_model_indices_C[0] = 0;
  for(int i=1; i<base_model_size; i++) base_model_indices_C[i]= INTEGER(base_model_indices_R)[i-1];
  data_store_constructor(X_in_C, y_in_C, weights_in_C, num_obs_in, num_var_in, base_model_indices_C, base_model_size);
  free(X_in_C); X_in_C = NULL;
  free(y_in_C); y_in_C = NULL;
  free(weights_in_C); weights_in_C = NULL;
  free(base_model_indices_C); base_model_indices_C = NULL;

  // setting model space prior
  void *vmax = vmaxget(); //just trying to get away from the R_alloc allocation for family_in
  const char * family_in = CHAR(STRING_ELT(msp_family_R,0));
  int params_length = Rf_length(msp_params_R);
  double * params_in = malloc(params_length*sizeof(double));
  for(int i=0; i<params_length; i++) params_in[i] = REAL(msp_params_R)[i];
  int trunc_in = INTEGER(trunc_R)[0];
  model_space_prior_constructor(family_in, params_in, params_length, trunc_in);
  vmaxset(vmax); //just trying to get away from the R_alloc allocation for family_in
  free(params_in); params_in=NULL;

  // log_bf_constructor(); //not yet

  //setting model fit
  model_fit_workspace_constructor();

  //setting hash key
  hash_key_parameters_constructor(data.bitrep_length);

  //setting rng
  // rng_constructor((unsigned long int) INTEGER(seed)[0]); //in mcmc_proposals

  //setting starting model
  int length = Rf_length(start_model_indices_R) + 1;
  int inc_coef = INTEGER(include_coef_R)[0];
  int inc_vcov = INTEGER(include_vcov_R)[0];
  int * ind = malloc(length*sizeof(int));
  ind[0] = 0;
  for(int i=1; i<length; i++) ind[i] = INTEGER(start_model_indices_R)[i-1];
  uint32_t * bitrep = int_to_uint32_t(ind, length, data.bitrep_length);
  model_t * start_model = model_bitrep_constructor(bitrep, data.bitrep_length, inc_coef, inc_vcov);
  div_t bucket_and_location;
  int needed_base_added = 0;
  for(int i=0; i<data.base_model_size; i++){
    bucket_and_location = div(data.base_model_indices[i], 32);
    if((start_model->bitrep[bucket_and_location.quot] & powers_of_2_uint32[bucket_and_location.rem]) == 0){
      start_model->bitrep[bucket_and_location.quot] |= powers_of_2_uint32[bucket_and_location.rem];
      start_model->size++;
      needed_base_added = 1;
    }
  }
  if(needed_base_added == 1) start_model->hash_key = hash_key_computer(start_model->bitrep, start_model->bitrep_length);
  free(ind); ind = NULL;
  free(bitrep); bitrep = NULL;

  // getting mcmc struct
  int burnin = INTEGER(n_burnin_R)[0];
  int thin = INTEGER(n_thin_R)[0];
  int draws = INTEGER(n_draws_R)[0];
  double prop_probs[3]; for(int i=0; i<3; i++) prop_probs[i] = REAL(proposal_probs_R)[i];
  mcmc_struct_constructor(start_model, burnin, thin, draws, prop_probs);
  model_destructor(&start_model);

  //getting hash table struct
  int ht_length = get_a_prime(INTEGER(hash_table_length_R)[0]);
  int ll_length = INTEGER(linked_list_length_R)[0];
  double mlf = REAL(max_load_factor_R)[0];
  hash_table = hash_table_constructor(ht_length, ll_length, mlf);

  //UNPROTECTING THE PROTECTED, the information is all on the C side now.
  UNPROTECT(17);

  // do the mcmc
  mcmc_all_draws(&hash_table);

  // release the memory that we can
  data_store_destructor();
  model_space_prior_destructor();
  // log_bf_destructor(); //not yet
  model_fit_workspace_destructor();
  hash_key_parameters_destructor();
  // rng_destructor();//in mcmc_proposals

  /* for now, this is what I want to output
   * models, size, rank, Rsq, log_prior, log_BF0, log_post_renormalization, mcmc_count, post_sampling, mcmc_id, mcmc_list
   */
  int out_elts_length = hash_table->hash_table_total_insertions;
  SEXP out = PROTECT(allocVector(VECSXP, 12));
  SEXP out_names = PROTECT(allocVector(STRSXP, 12));
  SEXP models = PROTECT(allocVector(VECSXP, out_elts_length));
  SEXP size = PROTECT(allocVector(INTSXP, out_elts_length));
  SEXP rank = PROTECT(allocVector(INTSXP, out_elts_length));
  SEXP mcmc_count = PROTECT(allocVector(INTSXP, out_elts_length));
  SEXP Rsq = PROTECT(allocVector(REALSXP, out_elts_length));
  SEXP residual_sd = PROTECT(allocVector(REALSXP, out_elts_length));
  SEXP log_prior = PROTECT(allocVector(REALSXP, out_elts_length));
  SEXP log_BF0 = PROTECT(allocVector(REALSXP, out_elts_length));
  SEXP log_post_renormalization = PROTECT(allocVector(REALSXP, out_elts_length));
  SEXP post_sampling = PROTECT(allocVector(REALSXP, out_elts_length));
  SEXP mcmc_id = PROTECT(allocVector(INTSXP, out_elts_length));
  SEXP mcmc_draws = PROTECT(allocVector(INTSXP, mcmc->num_draws));

  SEXP intarraySEXP;

  int mcmc_num_draws = mcmc->num_draws;

  for(int i=0; i<mcmc_num_draws; i++) INTEGER(mcmc_draws)[i] = mcmc->draws[i];
  mcmc_struct_destructor();
  int location = 0;
  int loc = 0;
  int * intarray = NULL;
  int intarray_length;
  for(int j=0; j<hash_table->hash_table_length; j++){
    for(int i=0; i<hash_table->linked_list_length; i++){
      if(hash_table->hash_table[loc+i].mod != NULL){
        INTEGER(size)[location] = hash_table->hash_table[loc+i].mod->size;
        INTEGER(rank)[location] = hash_table->hash_table[loc+i].mod->rank;
        INTEGER(mcmc_id)[location] = hash_table->hash_table[loc+i].mcmc_id;
        INTEGER(mcmc_count)[location] = hash_table->hash_table[loc+i].mcmc_count;
        REAL(residual_sd)[location] = hash_table->hash_table[loc+i].mod->residual_sd;
        REAL(Rsq)[location] = hash_table->hash_table[loc+i].mod->Rsq;
        REAL(log_prior)[location] = hash_table->hash_table[loc+i].mod->log_prior;
        REAL(log_BF0)[location] = hash_table->hash_table[loc+i].mod->log_BF0;
        REAL(log_post_renormalization)[location] = REAL(log_prior)[location] + REAL(log_BF0)[location];
        REAL(post_sampling)[location] = (double)(INTEGER(mcmc_count)[location])/(double)mcmc_num_draws;
        intarray = uint32_t_to_int(hash_table->hash_table[loc+i].mod->bitrep, hash_table->hash_table[loc+i].mod->bitrep_length, &intarray_length);
        intarraySEXP = PROTECT(allocVector(INTSXP, intarray_length));
        for(int k=0; k<intarray_length; k++) INTEGER(intarraySEXP)[k] = intarray[k];
        SET_VECTOR_ELT(models, location, intarraySEXP);
        UNPROTECT(1);
        free(intarray); intarray = NULL;
        model_destructor(&(hash_table->hash_table[loc+i].mod));
        location++;
      }else break;
    }
    loc += hash_table->linked_list_length;
  }
  free(hash_table); hash_table = NULL;

  double m = REAL(log_post_renormalization)[0];
  double s = 0;
  for(int i=1; i<out_elts_length; i++) m = fmax(m, REAL(log_post_renormalization)[i]);
  for(int i=0; i<out_elts_length; i++) s += exp(REAL(log_post_renormalization)[i]-m);
  s = log(s) + m;
  for(int i=0; i<out_elts_length; i++) REAL(log_post_renormalization)[i] -= s;

  SET_VECTOR_ELT(out, 0, models);
  SET_STRING_ELT(out_names, 0, mkChar("which"));

  SET_VECTOR_ELT(out, 1, log_BF0);
  SET_STRING_ELT(out_names, 1, mkChar("logmarg"));

  SET_VECTOR_ELT(out, 2, log_prior);
  SET_STRING_ELT(out_names, 2, mkChar("logprior"));

  SET_VECTOR_ELT(out, 3, log_post_renormalization);
  SET_STRING_ELT(out_names, 3, mkChar("logpost_RN"));

  SET_VECTOR_ELT(out, 4, post_sampling);
  SET_STRING_ELT(out_names, 4, mkChar("post_sampling"));

  SET_VECTOR_ELT(out, 5, size);
  SET_STRING_ELT(out_names, 5, mkChar("size"));

  SET_VECTOR_ELT(out, 6, rank);
  SET_STRING_ELT(out_names, 6, mkChar("rank"));

  SET_VECTOR_ELT(out, 7, Rsq);
  SET_STRING_ELT(out_names, 7, mkChar("Rsq"));

  SET_VECTOR_ELT(out, 8, residual_sd);
  SET_STRING_ELT(out_names, 8, mkChar("residual_sd"));

  SET_VECTOR_ELT(out, 9, mcmc_id);
  SET_STRING_ELT(out_names, 9, mkChar("mcmc_id"));

  SET_VECTOR_ELT(out, 10, mcmc_count);
  SET_STRING_ELT(out_names, 10, mkChar("mcmc_count"));

  SET_VECTOR_ELT(out, 11, mcmc_draws);
  SET_STRING_ELT(out_names, 11, mkChar("mcmc_draws"));

  Rf_setAttrib(out, R_NamesSymbol, out_names);

  // use the following when you want to pause before end and use leaks to check alloc/free
  // fprintf(stderr, "Press to end.\n");
  // getchar();

  PutRNGstate();

  UNPROTECT(14);
  return(out);

}

