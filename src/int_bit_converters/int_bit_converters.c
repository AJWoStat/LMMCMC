#include "int_bit_converters.h"

const uint32_t powers_of_2_uint32[32] = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288,1048576,2097152,4194304,8388608,16777216,33554432,67108864,134217728,268435456,536870912,1073741824,2147483648}; 

uint32_t * int_to_uint32_t(const int * x, const int x_length, const int out_length){
  uint32_t * out = (uint32_t *)calloc(out_length, sizeof(uint32_t));
  if(x==NULL || x_length==0) return(out);
  div_t bucket_and_location;
  for(int i= 0; i<x_length; i++){
    bucket_and_location = div(x[i], 32);
    out[bucket_and_location.quot] += powers_of_2_uint32[bucket_and_location.rem];
  }
  
  return(out);
}

int * uint32_t_to_int(const uint32_t * u, const int u_length, int * out_length){
  if(u==NULL || u_length==0) return((int*)NULL);
  int count = uint32_vec_count_bits(u, u_length);
  *out_length = count;
  if(count == 0) return((int*)NULL);
  int * out = malloc(count*sizeof(int));
  int out_location = 0;
  uint32_t value;
  for(int i=0; i<u_length; i++){
    value = u[i];
    for(int j=0; j<32; j++){
      if((value & powers_of_2_uint32[0]) > 0){
        out[out_location] = 32*i+j;
        out_location++;
      }
      value>>=1;
    }
  }
  return(out);
}

// int * uint32_t_to_int_with_tmp(const uint32_t * u, const int u_length, int * out_length, int ** tmp){
//   if(u==NULL || u_length==0) return((int*)NULL);
//   int count = 0;
//   int out_location = 0;
//   uint32_t value;
//   for(int i=0; i<u_length; i++){
//     value = u[i];
//     for(int j=0; j<32; j++){
//       if((value & powers_of_2_uint32[0]) > 0){
//         (*tmp)[out_location] = 32*i+j;
//         count++;
//         out_location++;
//       }
//       value>>=1;
//     }
//   }
//   
//   *out_length = count;
//   if(count==0) return((int*)NULL);
//   int * out = malloc(count*sizeof(int));
//   for(int i=0; i<count; i++) out[i] = (*tmp)[i];
//   
//   return(out);
// }

int uint32_vec_count_bits(const uint32_t *u, const int u_length){
  if(u==NULL || u_length==0) return(0);
  int count = 0;
  uint32_t value;
  for(int i=0; i<u_length; i++){
    value = u[i];
    for(int j=0; j<32; j++){
      if((value & powers_of_2_uint32[0]) > 0) count++;
      value>>=1;
    }
  }
  return(count);
}