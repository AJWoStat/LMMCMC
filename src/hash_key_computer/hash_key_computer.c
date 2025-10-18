#include"hash_key_computer.h"

//note use of rand() from stdlib to keep hash/rehash away from RNG from R

static struct hash_key_parameters_struct * hash_key_par = NULL;

void hash_key_parameters_constructor(const int vec_mult_length){
  hash_key_par = malloc(sizeof(struct hash_key_parameters_struct) + vec_mult_length*sizeof(uint64_t));
  hash_key_par->shift_right_xor[0] = 17;
  hash_key_par->shift_right_xor[1] = 11;
  hash_key_par->shift_right_xor[2] = 15;
  hash_key_par->shift_right_xor[3] = 14;
  hash_key_par->times_replace[0] = 0xed5ad4bb;
  hash_key_par->times_replace[1] = 0xac4c1b51;
  hash_key_par->times_replace[2] = 0x31848bab;
  
  int shift = (int)log2((double)(RAND_MAX)+1.00);
  int num_shift = 63/shift + 1;
  hash_key_par->vec_add = 0;
  for(int j=0; j<num_shift; j++){
    hash_key_par->vec_add += (((uint64_t)rand())<<(shift*j));
  }
  
  hash_key_par->vec_mult_length  = vec_mult_length;
  for(int i=0; i<vec_mult_length; i++){
    hash_key_par->vec_mult[i] = 0;
    for(int j=0; j<num_shift; j++){
      hash_key_par->vec_mult[i] += (((uint64_t)rand())<<(shift*j));
    }
  }
  return;
}

void hash_key_parameters_destructor(){
  if(hash_key_par!=NULL){free(hash_key_par); hash_key_par = NULL;}
}

uint32_t triple32inc(uint32_t u){
  /*
   * from https://github.com/skeeto/hash-prospector
   * by Christopher Wellons
   */
  u++;
  u ^= u >> hash_key_par->shift_right_xor[0];
  u *= hash_key_par->times_replace[0];
  u ^= u >> hash_key_par->shift_right_xor[1];
  u *= hash_key_par->times_replace[1];
  u ^= u >> hash_key_par->shift_right_xor[2];
  u *= hash_key_par->times_replace[2];
  u ^= u >> hash_key_par->shift_right_xor[3];
  return u;
}

uint32_t hash_key_computer(const uint32_t * u, const int u_length){
  /*
   * from https://arxiv.org/abs/1504.06804
   * by Mikkel Thorup
   */
  if(u_length != hash_key_par->vec_mult_length){
    Rprintf("Bad match, hashkey mutiplier parameter length does not match input bitrep length. Hash key 0 returned.");
    return(0);
  }
  uint64_t out = hash_key_par->vec_add;
  
  uint64_t tmp;
  for(int i=0; i<u_length; i++){
    tmp = (uint64_t)triple32inc(u[i]);
    out += hash_key_par->vec_mult[i]*tmp;
    
  }
  return((uint32_t)(out & 0xFFFFFFFF));
}


void hash_key_parameters_print(){
  Rprintf("Hash struct:\n Hash key xor: ");
  for(int i=0; i<4; i++) Rprintf("%u ", hash_key_par->shift_right_xor[i]);
  Rprintf("\nHash struct:\n Hash key times replace: ");
  for(int i=0; i<3; i++) Rprintf("%016llx ", hash_key_par->times_replace[i]);
  Rprintf("\nHash struct:\n Hash key add: %016llx\n", hash_key_par->vec_add);
  Rprintf("Hash struct:\n Hash key vec mult: ");
  for(int i=0; i<hash_key_par->vec_mult_length; i++) Rprintf("%016llx ", hash_key_par->vec_mult[i]);
  Rprintf("\n");
  
}
