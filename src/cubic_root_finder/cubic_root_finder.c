#include "cubic_root_finder.h"

//NOTE: RETURNS ALL REAL ROOTS, DOES NOT DEAL WITH COMPLEX ROOTS.
double * cubic_root_finder(double par[4], int * out_length){
  //get coefs
  double a = par[3]; //must be non-zero, but goes unchecked
  double b = par[2];
  double c = par[1];
  double d = par[0]; //must be non-zero, but goes unchecked
  //should really do some kind of a catch here
  //if a is 0, then go to a quadratic solver
  //if a is not 0 and d/a is 0, get 0 as a root and go to a quadratic solver
  //these would all have to account for the sizes of other coefficients
  //not going to code this up now, it should not be able to happen in the cases
  //that are needed for LMMCMC
  //but it might be worth either confirming or testing with edge cases of R^2
  
  
  //depressed cubic
  double root_adjust = b/(3.00*a);
  double p_1 = c/a;
  double p_2 = 3.00*root_adjust*root_adjust;
  double p = p_1-p_2;
  double q_1 = 2.00*root_adjust*root_adjust*root_adjust;
  double q_2 = root_adjust*p_1;
  double q_3 = d/a;
  double q = q_1-q_2+q_3;
  //discriminant of depressed cubic
  double delta_1 = 4.00*p*p*p;
  double delta_2 = 27.00*q*q;
  double delta = -(delta_1+delta_2);
  int out_depressed_length;
  double * out_depressed;
  //test number of real roots
  if(double_is_equal(p_1, p_2) && double_is_equal(q_2-q_1, q_3)){  
    //multiple roots scenario
    out_depressed_length = 3;
    out_depressed = malloc(out_depressed_length*sizeof(double));
    //triple root
    for(int i=0; i<3; i++) out_depressed[i] = 0.00;
  }else if(double_is_equal(delta_1, -delta_2)==1){
    //multiple roots scenario
    out_depressed_length = 3;
    out_depressed = malloc(out_depressed_length*sizeof(double));
    //one double root and one single root
    out_depressed[0] = 3.00*q/p;
    out_depressed[1] = -out_depressed[0]/2.00;
    out_depressed[2] = out_depressed[1];
  }else if(delta<0){
    //one distinct real root
    out_depressed_length = 1;
    out_depressed = malloc(out_depressed_length*sizeof(double));
    if(double_is_equal(p_1, p_2)==1){
      out_depressed[0] = -copysign(1.00, q)*pow(fabs(q),1.00/3.00);
    }else if(p>0){
      out_depressed[0] = -2.00*sqrt(p/3.00)*sinh(1.00/3.00*asinh((3.00*q)/(2.00*p)*sqrt(3.00/p)));
    }else{
      out_depressed[0] = -2.00*copysign(1.00, q)*sqrt(-p/3.00)*cosh(1/3.00*acosh((-3.00*fabs(q))/(2.00*p)*sqrt(-3.00/p)));
    }
  }else if(delta>0){
    //three distinct real roots
    out_depressed_length = 3;
    out_depressed = malloc(out_depressed_length*sizeof(double));
    for(int i=0; i<3; i++) out_depressed[i] = 2.00*sqrt(-p/3.00)*cos(1.00/3.00*acos((3.00*q)/(2.00*p)*sqrt(-3.00/p)) - 2.00*M_PI/3.00*i);
  }
  (*out_length) = out_depressed_length;
  double * out = malloc((*out_length)*sizeof(double));
  for(int i=0; i<(*out_length); i++) out[i] = out_depressed[i] - root_adjust;
  return(out);
}