#include "cubic_root_finder.h"

struct cubic_root_finder_struct crfs = {
  1.00, 1.00, 1.00, 1.00,
  'N', 'N',
  3, 
  0.00, 1.00, 0.00, 
  0.00, 0.00, 1.00, 
  1.00, 1.00, 1.00, 
  3,
  0.00, 0.00, 0.00, 
  0.00, 0.00, 0.00, 
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
  3, 3,
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 
  6,
  0
  };

void cubic_root_finder(double * par, double * output){
  /*
   * structure of the output is a length six double array 
   * the 0-2 elements are the real parts of the solutions
   * the 3-5 elements are the corresponding imaginary parts
   */
  
  //get coefs
  crfs.a = par[3]; //must be non-zero, but goes unchecked
  crfs.b = par[2];
  crfs.c = par[1];
  crfs.d = par[0]; 
  crfs.A[0] = 0.00; crfs.A[1] = 1.00; crfs.A[2] = 0.00;
  crfs.A[3] = 0.00; crfs.A[4] = 0.00; crfs.A[5] = 1.00;
  crfs.A[6] = -crfs.d/crfs.a; crfs.A[7] = -crfs.c/crfs.a; crfs.A[8] = -crfs.b/crfs.a;
  
  F77_NAME(dgeev)(&crfs.JOBVL, &crfs.JOBVR, 
               &crfs.N, crfs.A, &crfs.LDA,
               crfs.WR, crfs.WI, 
               crfs.VL, &crfs.LDVL, 
               crfs.VR, &crfs.LDVR,
               crfs.WORK, &crfs.LWORK,
               &crfs.INFO
                    FCONE FCONE
  );
  
  for(int i=0, j=3; i<3; i++, j++){
    output[i] = crfs.WR[i];
    output[j] = crfs.WI[j];
  }
  
  return;
}