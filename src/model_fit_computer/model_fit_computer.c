#include "model_fit_computer.h"

static struct model_fit_workspace_struct ws;

void model_fit_workspace_constructor(){
  //for general pre-allocation
  ws.XtX = (double*)calloc(data.num_var*data.num_var, sizeof(double));
  for(int i=0; i<data.num_var*data.num_var; i++) ws.XtX[i] = 0.0;
  ws.Xty = (double*)calloc(data.num_var, sizeof(double));
  for(int i=0; i<data.num_var; i++) ws.Xty[i] = 0.0;
  ws.XtX_piv = (double*)calloc(data.num_var*data.num_var, sizeof(double));
  for(int i=0; i<data.num_var*data.num_var; i++) ws.XtX_piv[i] = 0.0;
  ws.Xty_piv = (double*)calloc(data.num_var, sizeof(double));
  for(int i=0; i<data.num_var; i++) ws.Xty_piv[i] = 0.0;
  ws.model_indices = (int*)calloc(data.num_var, sizeof(int));
  for(int i=0; i<data.num_var; i++) ws.model_indices[i] = 0;

  //for using with LAPACK function dpstrf
  ws.uplo_lapack = 'U';
  ws.n_lapack = (La_INT) data.num_var;
  ws.a_lapack_ptr = ws.XtX;
  ws.lda_lapack = (La_INT) data.num_var;
  ws.piv_lapack_ptr = (La_INT*)calloc(data.num_var, sizeof(La_INT));;
  La_INT rank_lapack = 0;
  ws.work_lapack_ptr = (double*)calloc(2*data.num_var, sizeof(double));;
  ws.info_lapack = 0;
  char dlamch_e = 'e';
  ws.tol_lapack = F77_NAME(dlamch)(&dlamch_e FCONE)*data.sqrt_sum_sq_weights;



  //for using blas dtrsv
  ws.uplo_blas = 'U';
  ws.trans_blas = 'T';
  ws.diag_blas = 'N';
  ws.n_blas = (BLAS_INT)data.num_var;
  ws.a_blas_ptr = ws.XtX;
  ws.lda_blas = (BLAS_INT)data.num_var;
  ws.x_blas_ptr = ws.Xty_piv;
  ws.inc_x_blas = 1;

  //just for convenience
  ws.SST = 0;
  ws.SSE = data.base_model_SSE;

  //other stuff later

}
void model_fit_workspace_destructor(){
  free(ws.XtX); ws.XtX = NULL;
  free(ws.Xty); ws.Xty = NULL;
  free(ws.XtX_piv); ws.XtX_piv = NULL;
  free(ws.Xty_piv); ws.Xty_piv = NULL;
  free(ws.piv_lapack_ptr); ws.piv_lapack_ptr = NULL;
  free(ws.work_lapack_ptr); ws.work_lapack_ptr = NULL;
  free(ws.model_indices); ws.model_indices = NULL;
}

void model_fit_computer(model_t * mod){
  int model_indices_size;
  int * model_indices = uint32_t_to_int(mod->bitrep, mod->bitrep_length, &model_indices_size);
  for(int i=0; i<model_indices_size; i++) ws.model_indices[i] = model_indices[i];
  free(model_indices); model_indices = NULL;

  for(int j=0; j<mod->size; j++){
    for(int i=0; i<mod->size; i++){ //just taking upper triangle should be okay (it is all that is needed), but copying whole thing for error testing
      ws.XtX[i+j*mod->size] = data.XtX[ws.model_indices[i]+ws.model_indices[j]*data.num_var];
      // ws.XtX[j+i*mod->size] = ws.XtX[i+j*mod->size];
    }
    ws.Xty[j] = data.Xty[ws.model_indices[j]];
  }
  //replaces XtX_base with U where U^T*U = P^T*XtX_piv*P where P is pivoting matrix described by piv
  ws.n_lapack = (La_INT)model_indices_size;
  ws.lda_lapack = ws.n_lapack;
  F77_NAME(dpstrf)(&ws.uplo_lapack, &ws.n_lapack, &ws.a_lapack_ptr[0], &ws.lda_lapack, &ws.piv_lapack_ptr[0], &ws.rank_lapack, &ws.tol_lapack, ws.work_lapack_ptr, &ws.info_lapack FCONE);
  for(int i=0; i<(int)(ws.rank_lapack); i++) ws.Xty_piv[i] = ws.Xty[(int)(ws.piv_lapack_ptr[i])-1];
  // use dtrsv once just to get what we need to compute Rsq (no coef right now, add second call to get coef later)
  // replaces Xty_piv with (U^T)^(-1)*Xty_piv, recall XtX was replaced before with U by dpstrf
  ws.n_blas = (BLAS_INT) ws.rank_lapack;
  ws.lda_blas = ws.n_blas;
  F77_NAME(dtrsv)(&ws.uplo_blas, &ws.trans_blas, &ws.diag_blas, &ws.n_blas, &ws.a_blas_ptr[0], &ws.lda_blas, &ws.x_blas_ptr[0], &ws.inc_x_blas FCONE FCONE FCONE);
  ws.SST = F77_NAME(ddot)(&ws.n_blas, ws.x_blas_ptr, &ws.inc_x_blas, ws.x_blas_ptr, &ws.inc_x_blas);

  mod->Rsq = ws.SST/data.intercept_only_model_SSE; //if we want to use BF relative to base and not to intercept only, that is handled in the logBF function
  mod->residual_sd = sqrt((data.intercept_only_model_SSE-ws.SST)/((double)(data.num_obs - (int)ws.rank_lapack)));
  mod->rank = (int)ws.rank_lapack;

  
  mod->log_prior = model_space_prior.log_prior_models[mod->size - data.base_model_size];
  
  mod->log_BF0 = get_log_bf(mod);
  
  mod->shrinkage_factor = get_shrinkage_factor(mod, 1.00);
  
  
  //we need to finish with the coefs and such
  



}
