#include "data_store_struct.h"

struct data_store_struct data;

void data_store_constructor(double * X_in, double * y_in, double * weights_in, int num_obs_in, int num_var_in, int * base_model_indices_in, int base_model_size_in){
  data.num_obs = num_obs_in;
  data.num_var = num_var_in;
  data.base_model_size = base_model_size_in;
  data.num_test_var = data.num_var - data.base_model_size;
  data.base_model_size = base_model_size_in;
  data.XtX = (double*)calloc(data.num_var*data.num_var, sizeof(double));
  for(int i=0; i<data.num_var*data.num_var; i++) data.XtX[i] = 0;
  data.Xty = (double*)calloc(data.num_var, sizeof(double));
  for(int i=0; i<data.num_var; i++) data.Xty[i] = 0;
  data.weights = (double*)calloc(data.num_obs, sizeof(double));
  for(int i=0; i<data.num_obs; i++) data.weights[i] = weights_in[i];
  data.sqrt_sum_sq_weights = 0;
  for(int i=0; i<data.num_obs; i++) data.sqrt_sum_sq_weights += data.weights[i]*data.weights[i];
  data.sqrt_sum_sq_weights = sqrt(data.sqrt_sum_sq_weights);
  data.base_model_indices = (int*)calloc(data.base_model_size, sizeof(int));
  for(int i=0; i<data.base_model_size; i++) data.base_model_indices[i] = base_model_indices_in[i];
  data.bitrep_length = (data.num_var-1)/32 + 1;
  data.base_model_bitrep = int_to_uint32_t(data.base_model_indices, data.base_model_size, data.bitrep_length);
  data.test_var_indices = (int*)calloc(data.num_test_var, sizeof(int));
  int test_location = 0;
  int flag = 0;
  for(int i=0; i<data.num_var; i++){
    flag = 0;
    for(int j=0; j<data.base_model_size; j++){
      if(i == data.base_model_indices[j]){
        flag=1;
        break;
      }
    }
    if(flag==0){
      data.test_var_indices[test_location] = i;
      test_location++;
    }
  }

  // pre-processing data
  //multipy rows/observations by weights (unfortunately no hadamard product in blas)
  double * X_work = (double*)calloc(data.num_obs*data.num_var, sizeof(double));
  for(int j=0, k=0; j<data.num_var; j++, k+=data.num_obs){
    vec_hadamard_product(data.num_obs, &X_work[k], &X_in[k], data.weights);
  }

  double * y_work = (double*)calloc(data.num_obs, sizeof(double));
  vec_hadamard_product(data.num_obs, y_work, y_in, data.weights);

  //orthogonalize columns of X to first column of for numerical stability and make inner product or each column equal to sum_sq_weights
  BLAS_INT inc_blas = 1;
  BLAS_INT num_obs_blas = (BLAS_INT) data.num_obs;
  double sum_sq_weights = (data.sqrt_sum_sq_weights*data.sqrt_sum_sq_weights);
  double alpha_blas = data.sqrt_sum_sq_weights/F77_NAME(dnrm2)(&num_obs_blas, &X_work[0], &inc_blas);
  F77_NAME(dscal)(&num_obs_blas, &alpha_blas, &X_work[0], &inc_blas); //normalizing first column

  for(int j=1, k=data.num_obs; j<data.num_var; j++, k+=data.num_obs){
    alpha_blas = -F77_NAME(ddot)(&num_obs_blas, &X_work[0], &inc_blas, &X_work[k], &inc_blas)/sum_sq_weights;
    F77_NAME(daxpy)(&num_obs_blas, &alpha_blas, &X_work[0], &inc_blas, &X_work[k], &inc_blas);
    alpha_blas = data.sqrt_sum_sq_weights/F77_NAME(dnrm2)(&num_obs_blas, &X_work[k], &inc_blas);
    F77_NAME(dscal)(&num_obs_blas, &alpha_blas, &X_work[k], &inc_blas);
  }

  //same orthogonalizing and and normalizing for y - makes SSE for intercept only model into sum of squared weights
  alpha_blas = -F77_NAME(ddot)(&num_obs_blas, &X_work[0], &inc_blas, y_work, &inc_blas)/sum_sq_weights;
  F77_NAME(daxpy)(&num_obs_blas, &alpha_blas, &X_work[0], &inc_blas, y_work, &inc_blas);
  alpha_blas = data.sqrt_sum_sq_weights/F77_NAME(dnrm2)(&num_obs_blas, y_work, &inc_blas);
  F77_NAME(dscal)(&num_obs_blas, &alpha_blas, y_work, &inc_blas);



  //setting the intercept only SSE
  data.intercept_only_model_SSE = F77_NAME(ddot)(&num_obs_blas, y_work, &inc_blas, y_work, &inc_blas);

  //getting XtX and Xty
  char CblasUpper = 'U';
  char CblasTrans = 'T';
  BLAS_INT num_var_blas = (BLAS_INT) data.num_var;
  alpha_blas = 1.00;
  double beta_blas = 0.00;
  F77_NAME(dsyrk)(&CblasUpper, &CblasTrans, &num_var_blas, &num_obs_blas, &alpha_blas, X_work, &num_obs_blas, &beta_blas, data.XtX, &num_var_blas FCONE FCONE);
  F77_NAME(dgemv)(&CblasTrans, &num_obs_blas, &num_var_blas, &alpha_blas, X_work, &num_obs_blas, y_work, &inc_blas, &beta_blas, data.Xty, &inc_blas FCONE);

  fit_base_model();



  free(X_work); X_work = NULL;
  free(y_work); y_work = NULL;
  return;
}

void data_store_destructor(){
  free(data.base_model_indices); data.base_model_indices = NULL;
  free(data.test_var_indices); data.test_var_indices = NULL;
  free(data.base_model_bitrep); data.base_model_bitrep = NULL;
  free(data.XtX); data.XtX = NULL;
  free(data.Xty); data.Xty = NULL;
  free(data.weights); data.weights = NULL;
}

void vec_hadamard_product(const int n, double * out, const double * x, const double * y){
  for(int i=0; i<n; i++) out[i] = x[i]*y[i];
  return;
}

void fit_base_model(){ //returns rank of base model design matrix
  double * XtX_base = (double*)calloc(data.base_model_size*data.base_model_size, sizeof(double));
  double * Xty_base = (double*)calloc(data.base_model_size, sizeof(double));

  for(int j=0; j<data.base_model_size; j++){
    for(int i=0; i<data.base_model_size; i++){
      XtX_base[data.base_model_size*j+i] = data.XtX[data.num_var*data.base_model_indices[j] + data.base_model_indices[i]];
    }
    Xty_base[j] = data.Xty[data.base_model_indices[j]];
  }


  char uplo = 'U';
  char * uplo_ptr = &uplo;
  La_INT * piv = (La_INT*)calloc(data.base_model_size,sizeof(La_INT));
  double * work = malloc(2*data.base_model_size*sizeof(double));
  La_INT rank = 0;
  La_INT * rank_ptr = &rank;
  La_INT info = 0;
  La_INT * info_ptr = &info;
  La_INT size = (La_INT)data.base_model_size;
  La_INT * size_ptr = &size;
  double lapack_tol = 1e-12*data.sqrt_sum_sq_weights; // needs to be eps*sqrt(size of diag element of XtX after normalization), for us that is sqr(sum(weights^2))
  double * lapack_tol_ptr = &lapack_tol;

  //replaces XtX_base with U where U^T*U = P^T*XtX_piv*P where P is pivoting matrix described by piv
  F77_NAME(dpstrf)(uplo_ptr, size_ptr, XtX_base, size_ptr, piv, rank_ptr, lapack_tol_ptr, work, info_ptr FCONE);

  data.base_model_rank = (int)rank;

  double * Xty_piv = malloc(data.base_model_rank*sizeof(double));
  for(int i=0; i<data.base_model_rank; i++) Xty_piv[i] = Xty_base[(int)piv[i]-1];

  // use dtrsv once just to get what we need to compute Rsq (no coef right now, add second call to get coef later)
  // replaces Xty_piv with (XtX_base^T)^(-1)*Xty_piv, recall XtX_base was replaced before
  char CblasUpper = 'U';
  char CblasTrans = 'T';
  char CblasNonUnit = 'N';
  BLAS_INT base_model_rank_blas = (BLAS_INT) data.base_model_rank;
  BLAS_INT base_model_size_blas = (BLAS_INT) data.base_model_size;
  BLAS_INT inc_blas = 1;
  F77_NAME(dtrsv)(&CblasUpper, &CblasTrans, &CblasNonUnit, &base_model_rank_blas, XtX_base, &base_model_size_blas, Xty_piv, &inc_blas FCONE FCONE FCONE);

  double SST_base = F77_NAME(ddot)(&base_model_rank_blas , Xty_piv, &inc_blas, Xty_piv, &inc_blas);

  double SSE_base = data.intercept_only_model_SSE - SST_base;

  data.base_model_SSE = SSE_base;

  free(Xty_piv); Xty_piv = NULL;
  free(XtX_base); XtX_base = NULL;
  free(Xty_base); Xty_base = NULL;
  free(piv); piv= NULL;
  free(work); work= NULL;

  return;
}
