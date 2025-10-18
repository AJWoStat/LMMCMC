#ifndef MODEL_FIT_COMPUTER_H
#define MODEL_FIT_COMPUTER_H

#include <stdlib.h>
#include <math.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include "../model_struct/model_struct.h"
#include "../int_bit_converters/int_bit_converters.h"
#include "../data_store_struct/data_store_struct.h"
#include "../model_space_prior/model_space_prior.h"
#include "../log_bf_integrator/log_bf_integrator.h"


struct model_fit_workspace_struct{
  //to store in pre-allocated memory and reuse, initialized to have all 0s
  double * XtX;
  double * Xty;
  int * model_indices;
  double * XtX_piv; //this is just here if we need to get vcov, everything else can be done in XtX
  double * Xty_piv;

  //for using with LAPACK function dpstrf
  char uplo_lapack;
  La_INT n_lapack;
  double * a_lapack_ptr;
  La_INT lda_lapack;
  La_INT * piv_lapack_ptr;
  La_INT rank_lapack;
  double * work_lapack_ptr;
  La_INT info_lapack;
  double tol_lapack;

  //for using blas dtrsv
  char uplo_blas;
  char trans_blas;
  char diag_blas;
  BLAS_INT n_blas;
  double * a_blas_ptr;
  BLAS_INT lda_blas;
  double * x_blas_ptr;
  BLAS_INT inc_x_blas;

  //some values for convenience
  double SST;
  double SSE;


};

//static struct model_fit_workspace_struct ws;

void model_fit_workspace_constructor();
void model_fit_workspace_destructor();

void model_fit_computer(model_t * mod);

#endif
