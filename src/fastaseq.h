/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#ifndef _uvaia_fastaseq_h_
#define _uvaia_fastaseq_h_

#include <biomcmc.h>
#include "utils.h" 

typedef struct fastaseq_struct* fastaseq_t;
typedef struct cluster_struct* cluster_t;
typedef struct readfasta_struct* readfasta_t;

struct fastaseq_struct
{
  char **nn, *name, *seq;
  int n_nn;
  double score[2];
  size_t nchars;
};

struct cluster_struct
{
  fastaseq_t *fs;
  int n_fs, step;
  double mindist;
  size_t trim, nchars;
  char *reference;
};

struct readfasta_struct
{
  file_compress_t seqfile;
  char *line_read, *next_name;
  char *name, *seq; /*! \brief sequence and header */
  size_t linelength, seqlength;
  bool newseq;
};

int compare_fastaseq (const void *a, const void *b); /* \brief more neighours first, compare_score() if tie */
int compare_fastaseq_score (const void *a, const void *b);
fastaseq_t new_fastaseq (void);
void del_fastaseq (fastaseq_t fs);
void update_fasta_seq (fastaseq_t to, char **seq, char **name, size_t nchars, double *score);

cluster_t new_cluster (char *seq, size_t nchars, int mindist, size_t trim);
void del_cluster (cluster_t clus);
void check_seq_against_cluster (cluster_t clust, char **seq, char **name, size_t nchars);
void add_seq_to_cluster (cluster_t clust, int idx, char **seq, char **name, size_t nchars, double *score);
int merge_clusters (cluster_t clust1, cluster_t clust2); 
int compact_cluster (cluster_t clust);

void save_cluster_to_xz_file (cluster_t *clust, int n_clust, const char* filename);
void save_cluster_to_gz_file (cluster_t *clust, int n_clust, const char* filename);
void save_neighbours_to_xz_file (cluster_t *clust, int n_clust, const char* filename);
void save_neighbours_to_gz_file (cluster_t *clust, int n_clust, const char* filename);

readfasta_t new_readfasta (const char *seqfilename);
int readfasta_next (readfasta_t rfas);
void del_readfasta (readfasta_t rfas);

#endif
