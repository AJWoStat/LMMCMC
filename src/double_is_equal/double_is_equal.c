#include"double_is_equal.h"

double rel_tol = 1.00e-8;//pow(DBL_EPSILON, 0.50);
double abs_tol = 1.00e-16;//pow(DBL_EPSILON, 1.00);

int double_is_equal(const double a, const double b){
  if(a==0.0){
    if(fabs(b)>abs_tol) return(0);
  }else if(b==0.0){
    if(fabs(a)>abs_tol) return(0);
  }else{
    if(fabs(a-b)>(fmax(fabs(a),fabs(b))*rel_tol)) return(0);
  }
  return(1);
}


