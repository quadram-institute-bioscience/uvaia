/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#ifndef _uvaia_utils_h_
#define _uvaia_utils_h_

#include <biomcmc.h>
#include <gap_affine/affine_wavefront_align.h>
#include "kseq.h"


void add_reference_genome_to_char_vectors (char *name, char *s, unsigned l, char_vector cv_seq, char_vector cv_name);
double query_genome_against_char_vectors (char *name, char *s, unsigned l, char_vector cv_seq, char_vector cv_name, int nbest, int nmax, int **idx, int *n_idx);
void print_score_header (void);
void save_sequences (const char *filename, int *idx, int n_idx, char_vector seq, char_vector name);

#endif
