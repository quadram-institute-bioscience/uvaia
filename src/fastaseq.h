/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#ifndef _uvaia_fastaseq_h_
#define _uvaia_fastaseq_h_

#include <biomcmc.h>

typedef struct fastaseq_struct fastaseq_t;
typedef struct readfasta_struct readfasta_t;

struct fastaseq_struct
{
  char **nn, *name, *seq;
  int n_nn, score;
  size_t nchars;
}

struct readfasta_struct
{
  file_compress_t seqfile;
  char *line, *line_read, *next_name;
  size_t linelength;
  fastaseq_t fs;
}
#endif
