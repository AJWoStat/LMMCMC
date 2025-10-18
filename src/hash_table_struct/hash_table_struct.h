#ifndef MODEL_HASH_TABLE_STRUCT_H
#define MODEL_HASH_TABLE_STRUCT_H

#include<stdio.h>
#include<stdint.h>
#include<stdlib.h>
#include<tgmath.h>
#include"../model_struct/model_struct.h"
#include"../model_fit_computer/model_fit_computer.h"
#include"../int_bit_converters/int_bit_converters.h"
#include"../get_a_prime/get_a_prime.h"
#include"../hash_key_computer/hash_key_computer.h"
#include"R_ext/Print.h"


struct linked_list_element_struct{
  int mcmc_count;
  int mcmc_id;
  model_t * mod;
};

typedef struct linked_list_element_struct linked_list_element_t;

struct hash_table_struct{
  int hash_table_length;
  int linked_list_length;
  int hash_table_total_length;
  int hash_table_linked_lists_used;
  int hash_table_total_insertions;
  double load_factor;
  double load_factor_increase;
  double max_load_factor;//this is going to be compared to total_insertions/hash_table_total_length
  linked_list_element_t hash_table[];
};

typedef struct hash_table_struct hash_table_t;
typedef hash_table_t * hash_table_ptr_t;

extern hash_table_t * hash_table;


// modify constructor and destructor
// need insertion operator, maybe a rehasher/resizer
// but the rehasher will have to rehash all of the models and remake all of the linked lists
// so this only matters if I am not resizing the linked lists
// do I want to resize the linked lists or do I want to rehash?
// the plan is to make the hash table twice the desired size (pointers are cheap)
// if I am using a good hash function there should be few collisions, so the linked lists should all be short
hash_table_t *  hash_table_constructor(const int hash_table_length_in, const int linked_list_length_in, const double max_load_factor_in);
void hash_table_destructor(hash_table_t ** ht);
int hash_table_insertion(hash_table_t ** ht, model_t * mod, int * hash_table_location_out, int * linked_list_location_out);
int hash_table_resize(hash_table_t ** ht, int length_new);
void hash_table_rehash(hash_table_t * ht, const int bitrep_length);
int hash_table_resize_insertion(hash_table_t ** ht_new, linked_list_element_t ll_elt);
void hash_table_resize_fail_reinsertion(hash_table_t ** ht, linked_list_element_t  ll_elt);

#endif
