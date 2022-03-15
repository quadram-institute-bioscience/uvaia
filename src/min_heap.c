/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "min_heap.h"

int compare_ratchet_element_score (int *a, int *b);

ratchet_t
new_ratchet (int n)
{
  int i,j;
  ratchet_t r = (ratchet_t) biomcmc_malloc (sizeof (struct ratchet_struct));
  r->sorted = false;
  r->n = n;
  r->current = 0;
  r->seq = (ratchet_element*) biomcmc_malloc (n * sizeof (ratchet_element));
  for (i = 0; i < n; i++) {
    r->seq[i].name = NULL;
    for (j = 0; j < 4; j++) r->seq[i].score[j] = 0;
  }
  return r;
}

void
del_ratchet (ratchet_t r)
{
  int i;
  if (!r) return;
  if (r->seq) {
    for (i = r->n - 1; i >= 0; i --) if (r->seq[i].name) free (r->seq[i].name);
    free (r->seq);
  }
  free (r);
}

int
compare_ratchet_element (const void *a, const void *b)
{
  return compare_ratchet_element_score(((ratchet_element*)a)->score, ((ratchet_element*)b)->score);
}

int
compare_ratchet_element_score (int *a, int *b)
{ // ratchet is worse -> to -> best score (so if b should go first, it should return positive since qsort is ascending order)
  int res = 0;
  for (int i = 0; (i < 4) && (!res); i++) res = a[i] - b[i]; // increasing, thus higher scores go after (should be #matches, #valid)
  return res;
}

/**** min heap from amburana ****/

heap_hash64 
new_heap_hash64 (int heap_size)
{
  int i;
  heap_hash64 pq = (heap_hash64) biomcmc_malloc (sizeof (struct heap64_struct)); 
  pq->n = 0;
  pq->heap_size = (heap_size < 2 ? 2: heap_size);
  pq->item = (hpq_item*) biomcmc_malloc ((pq->heap_size + 1) * sizeof (hpq_item));
  for (i = 0; i < pq->heap_size + 1; i++) { pq->item[i].freq = pq->item[i].id0 = pq->item[i].id = 0; pq->item[i].hash = 0ULL; }
  return pq;
}

void
del_heap_hash64 (heap_hash64 pq)
{
  if (!pq) return;
  if (pq->item) free (pq->item); 
  free (pq); 
}

hpq_item
heap_hash64_get_maximum (heap_hash64 pq) { return pq->item[1]; }

hpq_item
heap_hash64_remove_maximum (heap_hash64 pq) // a.k.a. pop() 
{
  hpq_item max = pq->item[1]; //the root is the maximum element
  if (pq->n == 0) return (hpq_item) {0,0,0,0};
  pq->item[1] = pq->item[pq->n]; // replace the element at the top with last element in the heap
  pq->n--; 
  heap_hash64_bubble_down (pq, 1);
  return max;
}

void 
heap_hash64_insert (heap_hash64 pq, hpq_item item) // struct, not a pointer
{
  if (heap_hash64_item_already_here (pq, item, 1)) return; // if exists, then increase count  
  if (pq->n == pq->heap_size) { // replace if smaller than maximum
    if (item.hash < pq->item[1].hash) { 
      pq->item[1].id0  = item.id0; // id0 never changes, but id is updated with most recent location
      pq->item[1].id   = item.id; 
      pq->item[1].freq = item.freq;
      pq->item[1].hash = item.hash;
      heap_hash64_bubble_down (pq, 1); 
    }
  }
  else { // regular insert (a.k.a. push)
    pq->n++;
    pq->item[pq->n].id0  = item.id0; 
    pq->item[pq->n].id   = item.id;
    pq->item[pq->n].freq = item.freq; // should be one
    pq->item[pq->n].hash = item.hash;
    heap_hash64_bubble_up (pq, pq->n);
  }
}

static bool
heap_hash64_item_already_here (heap_hash64 pq, hpq_item item, int p)
{ // must traverse both children since MaxHeap is not ordered like BST
  return false; // note from 2022.03.15: unfinished idea is to increase count (histsketch) but currently we discard this info (if item is already here)
  if (p > pq->n) return false; // it can be equal, since we have one extra element (item[1...n] )
  if (item.hash == pq->item[p].hash) { pq->item[p].id = item.id; pq->item[p].freq += item.freq; return true; }
  if (item.hash > pq->item[p].hash) return false; /* above, notice that only id is updated (id0 is leftmost location) */
  bool below = heap_hash64_item_already_here (pq, item, 2 * p); // 'left' child 
  if (!below) below = heap_hash64_item_already_here (pq, item, 2 * p + 1); // traverse right child if not found yet
  return below;
}

static void
heap_hash64_bubble_down (heap_hash64 pq, int p)
{//bubbles down an element into it's proper place in the heap
	int i, c = 2*p, max_index = p; //c= index of younger child, max_index = index of maximum child
  for (i = 0; i < 2; i++) if(c + i <= pq->n) {
    //check to see if the data at min_index is smaller than the data at the child
    if (pq->item[max_index].hash < pq->item[c+i].hash)	max_index = c + i;
  }
	if (max_index != p) {
    hpq_item tmp = pq->item[p];
    pq->item[p] = pq->item[max_index];
    pq->item[max_index] = tmp;
    heap_hash64_bubble_down (pq, max_index);
  }
}

static void 
heap_hash64_bubble_up (heap_hash64 pq, int index)
{ //bubbles up the last element of the heap to maintain heap structure
  int parent = (index == 1 ? -1 : (int)(index/2));
  if (parent == -1) return;	//if we are at the root of the heap, no parent
  //if the parent node has a smaller data value, we need to bubble up
  if (pq->item[parent].hash < pq->item[index].hash) {
    hpq_item tmp = pq->item[parent];
    pq->item[parent] = pq->item[index];
    pq->item[index] = tmp;
    heap_hash64_bubble_up (pq, parent);
  }
}

void
heap_hash64_finalise_heap_qsort (heap_hash64 pq)
{
  hpq_item tmp = pq->item[0];  pq->item[0] = pq->item[pq->n];  pq->item[pq->n] = tmp;   // heap goes from [1...n] and we want [0...n-1]
  qsort (pq->item, pq->n, sizeof (hpq_item), compare_hpq_item_increasing);
  if (pq->n < pq->heap_size - 1) { 
    pq->item = (hpq_item*) biomcmc_realloc ((hpq_item*) pq->item, pq->n * sizeof (hpq_item));
    pq->heap_size = pq->n; 
  }
}

static int
compare_hpq_item_increasing (const void *a, const void *b)
{
  if ( ((hpq_item *)b)->hash < ((hpq_item *)a)->hash ) return 1;
  if ( ((hpq_item *)b)->hash > ((hpq_item *)a)->hash ) return -1;
  return 0;
}
