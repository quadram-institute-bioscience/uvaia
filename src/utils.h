/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#ifndef _uvaia_utils_h_
#define _uvaia_utils_h_

#include <biomcmc.h>
#include <gap_affine/affine_wavefront_align.h>
#include "kseq.h"

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

void uvaia_keep_only_valid_sequences (alignment aln, double ambiguity, bool check_aligned);
double query_genome_against_char_vectors (char *name, char *s, unsigned l, char_vector cv_seq, char_vector cv_name, 
                                          int nbest, int nmax, int **idx, int *n_idx, size_t trim);
void print_score_header (void);
void save_sequences (const char *filename, int *idx, int n_idx, char_vector seq, char_vector name);
char *return_query_aligned (int pattern_length, char* text, int text_length, edit_cigar_t* edit_cigar, mm_allocator_t* mm_allocator);
void upper_kseq (char *s, unsigned l);

void initialise_acgt (void);
bool is_site_acgt_distinct_pair (char s1, char s2); // relies on external call to initialise_acgt()
bool is_site_acgt (char s1); // relies on external call to initialise_acgt()

#endif
