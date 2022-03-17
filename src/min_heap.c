/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "min_heap.h"

#define N_SCORES 6  // don't forget to change the "return zero" in remove_worse() 

int compare_q_item_score (int *a, int *b);
static void heap_bubble_down (heap_t pq, int p);
static void  heap_bubble_up (heap_t pq, int index);

ratchet_t
new_ratchet (int n)
{
  int i,j;
  ratchet_t r = (ratchet_t) biomcmc_malloc (sizeof (struct ratchet_struct));
  r->sorted = false;
  r->n = n;
  r->current = 0;
  r->seq = (q_item*) biomcmc_malloc (n * sizeof (q_item));
  for (i = 0; i < n; i++) { r->seq[i].name = NULL; for (j = 0; j < N_SCORES; j++) {r->seq[i].score[j] = 0; } }
  return r;
}

void
del_ratchet (ratchet_t r)
{
  int i;
  if (!r) return;
  if (r->seq) { for (i = r->n - 1; i >= 0; i --) { if (r->seq[i].name) free (r->seq[i].name);} free (r->seq);  }
  free (r);
}

int
compare_q_item (const void *a, const void *b)
{
  return compare_q_item_score(((q_item*)a)->score, ((q_item*)b)->score);
}

int
compare_q_item_score (int *a, int *b)
{ // ratchet is worse -> to -> best score (so if b should go first, it should return positive since qsort is ascending order)
  int res = 0;
  for (int i = 0; (i < N_SCORES) && (!res); i++) res = b[i] - a[i]; // decreasing, thus higher scores go first (should be #matches, #valid)
  return res;
} // remember that min_heap keep "smallest" and border (root) is largest, which keeps being replaced (like "mismatches"); here we insert if returns +1

/**** min heap from amburana ****/

heap_t
new_heap_t (int heap_size)
{
  int i,j;
  heap_t pq = (heap_t) biomcmc_malloc (sizeof (struct heap_struct)); // some ppl suggested to avoid the type and sizeof(variable)
  pq->n = 0;
  pq->max_incompatible = 0xffffff; // starts with large value
  pq->heap_size = (heap_size < 2 ? 2: heap_size);
  pq->seq = (q_item*) biomcmc_malloc ((pq->heap_size + 1) * sizeof (q_item));
  for (i = 0; i < pq->heap_size + 1; i++) { 
    pq->seq[i].name = NULL; 
    for (j = 0; j < N_SCORES; j++) pq->seq[i].score[j] = 0;
  }
  return pq;
}

void
del_heap_t (heap_t pq)
{
  int i;
  if (!pq) return;
  if (pq->seq) {
    for (i = pq->heap_size - 1; i >= 0; i --) if (pq->seq[i].name) free (pq->seq[i].name);
    free (pq->seq);
  }
  free (pq); 
}

q_item
heap_get_worse (heap_t pq) { return pq->seq[1]; } // in a min_heap this is the largest amongst the lowest

q_item
heap_remove_worse (heap_t pq) // a.k.a. pop() maximum in a min_heap
{
  q_item max = pq->seq[1]; //the root is the maximum element (max mismatches, using our "min" logic)
  if (pq->n == 0) return (q_item) {.name=NULL, .score={0,0,0,0,0,0}};
  pq->seq[1] = pq->seq[pq->n]; // replace the element at the top with last element in the heap
  pq->n--; 
  heap_bubble_down (pq, 1);
  return max;
}

bool
heap_insert (heap_t pq, q_item item) // item is a struct, not a pointer
{ 
  if (pq->n == pq->heap_size) { // replace if smaller than maximum
    if (compare_q_item_score (item.score, pq->seq[1].score) < 0) { // item has better score than the worst in heap
      //pq->seq[1].name = (char*) biomcmc_realloc ((char*) pq->seq[1].name, strlen(item.name) * sizeof (char));
      //strcpy(pq->seq[1].name, item.name); // copy is better than relying on ptrs
      if (pq->seq[1].name) free (pq->seq[1].name);
      pq->seq[1].name = strdup (item.name);
      for (int j = 0; j < N_SCORES; j++) pq->seq[1].score[j] = item.score[j];
      heap_bubble_down (pq, 1); 
      return true;
    }
    else return false; // new item is not better than existing 
  }
  else { // regular insert (a.k.a. push)
    pq->n++;
    //pq->seq[pq->n].name = (char*) biomcmc_malloc (strlen(item.name) * sizeof (char));
    //strcpy (pq->seq[pq->n].name, item.name);
    pq->seq[pq->n].name = strdup (item.name);
    for (int j = 0; j < N_SCORES; j++) pq->seq[pq->n].score[j] = item.score[j];
    heap_bubble_up (pq, pq->n);
    return true; // new item was added
  }
}

static void
heap_bubble_down (heap_t pq, int p)
{//bubbles down an element into its proper place in the heap
	int i, c = 2*p, max_index = p; //c= index of younger child, max_index = index of maximum child
  for (i = 0; i < 2; i++) if(c + i <= pq->n) {
    //check to see if the data at min_index is smaller than the data at the child
    if (compare_q_item_score (pq->seq[max_index].score, pq->seq[c+i].score) < 0) max_index = c + i; // a has better score than b
  }
	if (max_index != p) {
    q_item tmp = pq->seq[p];
    pq->seq[p] = pq->seq[max_index];
    pq->seq[max_index] = tmp;
    heap_bubble_down (pq, max_index);
  }
}

static void 
heap_bubble_up (heap_t pq, int index)
{ //bubbles up the last element of the heap to maintain heap structure
  int parent = (index == 1 ? -1 : (int)(index/2));
  if (parent == -1) return;	//if we are at the root of the heap, no parent
  //if the parent node has a smaller data value, we need to bubble up
  if (compare_q_item_score (pq->seq[parent].score, pq->seq[index].score) < 0) { // parent has better score than index
    q_item tmp = pq->seq[parent];
    pq->seq[parent] = pq->seq[index];
    pq->seq[index] = tmp;
    heap_bubble_up (pq, parent);
  }
}

void
heap_finalise_heap_qsort (heap_t pq)
{
  q_item tmp = pq->seq[0];  pq->seq[0] = pq->seq[pq->n]; pq->seq[pq->n] = tmp; // heap goes from [1...n] and we want [0...n-1]
  qsort (pq->seq, pq->n, sizeof (q_item), compare_q_item);
  if (pq->n < pq->heap_size - 1) { // we can assume that n > heap_size elements are NULL, never used 
    pq->seq = (q_item*) biomcmc_realloc ((q_item*) pq->seq, pq->n * sizeof (q_item));
    pq->heap_size = pq->n; 
  }
}

