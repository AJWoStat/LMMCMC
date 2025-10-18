#include"hash_table_struct.h"

hash_table_t * hash_table;

hash_table_t * hash_table_constructor(const int hash_table_length_in, const int linked_list_length_in, const double max_load_factor_in){
  hash_table_t * ht;
  ht = malloc(sizeof(hash_table_t) + hash_table_length_in*linked_list_length_in*sizeof(linked_list_element_t));
  if(ht==NULL) fprintf(stderr,"hash table allocation failed\n");
  ht->hash_table_length = hash_table_length_in;
  ht->linked_list_length = linked_list_length_in;
  ht->hash_table_total_length = (ht->hash_table_length)*(ht->linked_list_length);
  ht->hash_table_linked_lists_used = 0;
  ht->hash_table_total_insertions = 0;
  ht->load_factor = 0.00;
  ht->load_factor_increase = 1.00/((double)ht->hash_table_total_length);
  ht->max_load_factor = max_load_factor_in;
  for(int i=0; i<ht->hash_table_total_length; i++){
    ht->hash_table[i].mcmc_count = 0;
    ht->hash_table[i].mcmc_id = -1;
    ht->hash_table[i].mod = NULL;
  }
  return(ht);
}

void hash_table_destructor(hash_table_t ** ht){
  int loc = 0;
  if((*ht)!=NULL){
    for(int j=0; j<(*ht)->hash_table_total_length; j++){
      for(int i=0; i<(*ht)->hash_table_total_length; i++){
        if((*ht)->hash_table[loc+i].mod != NULL){
          model_destructor(&((*ht)->hash_table[loc+i].mod));
        }
      }
      loc += ((*ht)->hash_table_total_length);
    }
    free((*ht)); (*ht) = NULL;
  }
  return;
}

int hash_table_insertion(hash_table_t ** ht, model_t * mod, int * hash_table_location_out, int * linked_list_location_out){
  //the first input is needed in case a resize has to happen and we need to replace the hash table
  //returns 0 if no resize needed and 1 if resize was needed
  //side effect here of changing the hash table and linked list location outputs in the arguments list, those inputs are not used and are only replaced
  //side effect of fitting or copying fit values for mod
  uint32_t htu32 = (mod->hash_key%((uint32_t)(*ht)->hash_table_length));
  int hash_table_location = (int)(mod->hash_key%((uint32_t)(*ht)->hash_table_length));
  int linked_list_location;
  int loc = hash_table_location*((*ht)->linked_list_length);
  for(int i=0; i<(*ht)->linked_list_length; i++){
    if((*ht)->hash_table[loc+i].mod == NULL){
      if(i==0) (*ht)->hash_table_linked_lists_used++; 
      (*ht)->hash_table_total_insertions++;
      (*ht)->load_factor += (*ht)->load_factor_increase;
      linked_list_location = i;
      model_fit_computer(mod);
      (*ht)->hash_table[loc+i].mod = model_copy_constructor(mod);
      break;
    }else if(model_is_equal((*ht)->hash_table[loc+i].mod, mod, 1)==1){
      linked_list_location = i;
      model_copy(mod, (*ht)->hash_table[loc+i].mod);
      break;
    }
  }
  int resize_flag=0;
  int new_length;
  int resize_fail;
  int num_resize;
  if(linked_list_location == ((*ht)->linked_list_length-1) || (*ht)->load_factor>=(*ht)->max_load_factor){
    resize_flag = 1;
    new_length = get_a_prime((*ht)->hash_table_length*2);
    resize_fail = hash_table_resize(ht, new_length);
    num_resize=1;
    while(resize_fail==1 && num_resize<5){
      new_length = get_a_prime(new_length+1);
      resize_fail = hash_table_resize(ht, new_length);
      num_resize++;
      if(resize_fail==1 && num_resize==5){
        hash_table_rehash(*ht, mod->bitrep_length);
        resize_fail = hash_table_resize(ht, new_length);
        num_resize=1;
      }
    }
  }
  *hash_table_location_out = hash_table_location;
  *linked_list_location_out = linked_list_location;
  return(resize_flag);
}

int hash_table_resize(hash_table_t ** ht, int length_new){
  hash_table_t * ht_new = hash_table_constructor(length_new, (*ht)->linked_list_length, (*ht)->max_load_factor);
  int resize_fail = 0;
  int loc = 0;
  for(int j=0; j<(*ht)->hash_table_length; j++){
    for(int i=0; i<(*ht)->linked_list_length; i++){
      if((*ht)->hash_table[loc+i].mod != NULL){
        resize_fail = hash_table_resize_insertion(&ht_new, (*ht)->hash_table[loc+i]);
        model_destructor(&((*ht)->hash_table[loc+i].mod));
        if(resize_fail==1){
          loc = 0;
          for(int k=0; k<ht_new->hash_table_length; k++){
            for(int l=0; l<ht_new->linked_list_length; l++){
              if(ht_new->hash_table[loc+l].mod != NULL){
                hash_table_resize_fail_reinsertion(ht, ht_new->hash_table[loc + l]);
                model_destructor(&(ht_new->hash_table[loc + l].mod));
              }else break;
            }
            loc += (ht_new->linked_list_length);
          }
          free(ht_new); ht_new = NULL;
          return(1);
        }
      }else break;
    }
    loc += ((*ht)->linked_list_length);
  }
  free(*ht); (*ht) = NULL;
  (*ht) = ht_new;
  return(0);
}

void hash_table_rehash(hash_table_t * ht, const int bitrep_length){
  hash_key_parameters_destructor();
  hash_key_parameters_constructor(bitrep_length);
  int loc = 0;
  for(int j=0; j<ht->hash_table_length; j++){
    for(int i=0; i<ht->linked_list_length; i++){
      if(ht->hash_table[loc+i].mod!=NULL){
        ht->hash_table[loc+i].mod->hash_key = hash_key_computer(ht->hash_table[loc+i].mod->bitrep, ht->hash_table[loc+i].mod->bitrep_length);  
      }else break;
    }
    loc += (ht->linked_list_length);
  }
  return;
}

int hash_table_resize_insertion(hash_table_t ** ht_new, linked_list_element_t  ll_elt){
  //returns 1 if insertion triggers a resize, 0 otherwise
  int hash_table_location = (int)(ll_elt.mod->hash_key%((uint32_t)(*ht_new)->hash_table_length));
  int htl = hash_table_location*((*ht_new)->linked_list_length);
  int linked_list_location;
  for(int i=0; i<(*ht_new)->linked_list_length; i++){
    if((*ht_new)->hash_table[htl+i].mod == NULL){
      if(i==0) (*ht_new)->hash_table_linked_lists_used++;
      (*ht_new)->hash_table_total_insertions++;
      (*ht_new)->load_factor += (*ht_new)->load_factor_increase;
      linked_list_location = i;
      (*ht_new)->hash_table[htl+i].mcmc_count = ll_elt.mcmc_count;
      (*ht_new)->hash_table[htl+i].mcmc_id = ll_elt.mcmc_id;
      (*ht_new)->hash_table[htl+i].mod = model_copy_constructor(ll_elt.mod);
      break;
    }
  }
  int resize_fail=0;
  if(linked_list_location == ((*ht_new)->linked_list_length-1) || (*ht_new)->load_factor>=(*ht_new)->max_load_factor){
    resize_fail = 1;
  }
  return(resize_fail);
}

void hash_table_resize_fail_reinsertion(hash_table_t ** ht, linked_list_element_t  ll_elt){
  //returns 1 if insertion triggers a resize, 0 otherwise
  int hash_table_location = (int)(ll_elt.mod->hash_key%((uint32_t)(*ht)->hash_table_length));
  int htl = hash_table_location*((*ht)->linked_list_length);
  for(int i=0; i<(*ht)->linked_list_length; i++){
    if((*ht)->hash_table[htl+i].mod == NULL){
      (*ht)->hash_table[htl+i].mcmc_count = ll_elt.mcmc_count;
      (*ht)->hash_table[htl+i].mcmc_id = ll_elt.mcmc_id;
      (*ht)->hash_table[htl+i].mod = model_copy_constructor(ll_elt.mod);
      break;
    }
  }
}
