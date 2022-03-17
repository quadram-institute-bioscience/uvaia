/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#ifndef _uvaia_min_heap_h_
#define _uvaia_min_heap_h_

#include <biomcmc.h>
#include "utils.h"

typedef struct ratchet_struct* ratchet_t;
typedef struct heap_struct* heap_t;

typedef struct
{
  int score[6];
  char *name;
} q_item;

struct ratchet_struct
{
  int n, current;
  q_item *seq;
  bool sorted;
};

struct heap_struct {
  q_item *seq; 
  int heap_size, n;
  int max_incompatible; // unrelated to min_heap but used to truncate dist calculation
};

int compare_q_item_score (int *a, int *b);
heap_t new_heap_t (int heap_size);
void del_heap_t (heap_t pq);
q_item heap_get_worse (heap_t pq); // in a min_heap this is the largest amongst the lowest
q_item heap_remove_worse (heap_t pq); // a.k.a. pop() maximum in a min_heap
bool heap_insert (heap_t pq, q_item item); // item is a struct, not a pointer
void heap_finalise_heap_qsort (heap_t pq);

#endif
