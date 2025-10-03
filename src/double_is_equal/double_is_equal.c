#include"double_is_equal.h"

int double_is_equal(const double a, const double b){
  double rel_tol = pow(DBL_EPSILON, 0.85);
  double abs_tol = fmax(pow(DBL_EPSILON, 1.25), DBL_MIN);
  if(a==0.0){
    if(fabs(b)>abs_tol) return(0);
  }else if(b==0.0){
    if(fabs(a)>abs_tol) return(0);
  }else{
    if(fabs(a-b)>(fmax(fabs(a),fabs(b))*rel_tol)) return(0);
  }
  return(1);
}


