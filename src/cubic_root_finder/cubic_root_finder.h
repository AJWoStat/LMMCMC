#ifndef CUBIC_ROOT_FINDER_H
#define CUBIC_ROOT_FINDER_H

#include<stdlib.h>
#include<math.h>
#include"R_ext/Lapack.h"


struct cubic_root_finder_struct{
  double a;
  double b;
  double c;
  double d;
  char JOBVL, JOBVR;
  La_INT N;
  double A[9];
  La_INT LDA;
  double WR[3], WI[3];
  double VL[9], VR[9];
  La_INT LDVL, LDVR;
  double WORK[18];
  La_INT LWORK;
  La_INT INFO;
};

void cubic_root_finder(double * par, double * output);
/*
 * structure of the output is a length six double array 
 * the 0-2 elements are the real parts of the solutions
 * the 3-5 elements are the corresponding imaginary parts
 */

#endif