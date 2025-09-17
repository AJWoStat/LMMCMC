#include"double_is_equal.h"

static const double tol_die = 1e-10;

int double_is_equal(const double a, const double b){
  if((2.00*fabs(a-b)/(fabs(a)+fabs(b)))>tol_die) return(0);
  return(1);
}


