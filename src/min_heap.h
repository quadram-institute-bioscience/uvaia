/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#ifndef _uvaia_min_heap_h_
#define _uvaia_min_heap_h_

#include <biomcmc.h>
#include "utils.h"

typedef struct ratchet_struct* ratchet_t;

typedef struct
{
  int score[4];
  char *name;
} ratchet_element;

struct ratchet_struct
{
  int n, current;
  ratchet_element *seq;
  bool sorted;
};

#endif
