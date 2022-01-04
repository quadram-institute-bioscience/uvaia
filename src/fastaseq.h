/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#ifndef _uvaia_fastaseq_h_
#define _uvaia_fastaseq_h_

#include <biomcmc.h>

typedef struct fastaseq_struct* fastaseq_t;
typedef struct cluster_struct* cluster_t;
typedef struct readfasta_struct* readfasta_t;

struct fastaseq_struct
{
  char **nn, *name, *seq;
  int n_nn, score;
  size_t nchars;
};

struct cluster_struct
{
  fastaseq_t fs;
  int n_fs;
};

struct readfasta_struct
{
  file_compress_t seqfile;
  char *line_read, *next_name;
  char *name, *seq; /*! \brief sequence and header */
  size_t linelength, seqlength;
  bool newseq;
};

fastaseq_t new_fastaseq (void);
void del_fastaseq (fastaseq_t fs);
void update_fasta_seq (fastaseq_t to, char **seq, char **name, size_t nchars, int score);

readfasta_t new_readfasta (const char *seqfilename);
int readfasta_next (readfasta_t rfas);
void del_readfasta (readfasta_t rfas);

#endif
