#ifndef LOG_BF_COMPUTER_H
#define LOG_BF_COMPUTER_H

#include<stdlib.h>
#include<math.h>
#include"../model_struct/model_struct.h"
#include"../data_store_struct/data_store_struct.h"
// #include"model_struct.h"
// #include"data_store_struct.h"

//lots of work needs to be done here for getting the integrals, stealing code form SearchLM is a good idea
//for now, it is just a simple function that computes the standard g-over-n prior with g=1

// struct log_bf_workspace_struct{
//   int scale_by_effective_number_of_parameters_flag; //1 if you do, 0 if you don't
//   char * type;
//   int parameters_length;
//   double * parameters;
// };


double log_bf_computer(const model_t * mod);
double coef_shrinkage_factor_computer(const model_t * mod);
double coef_vcov_shrinkage_factor_computer(const model_t * mod);

#endif
