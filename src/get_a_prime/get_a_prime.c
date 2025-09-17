#include "get_a_prime.h"

int is_prime(int n){
  int rem = 1;
  int div = 2;
  while(rem>0 && div*div<=n){
    rem = n%div;
    div++;
  }
  if(rem>0) rem=1;
  return(rem);
}

int get_a_prime(int lower_bound){
  int out = lower_bound+1;
  if(out%2 == 0) out++;
  while(!is_prime(out)) out+=2;
  return(out);
}


