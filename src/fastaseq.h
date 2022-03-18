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
typedef struct query_struct* query_t;

struct fastaseq_struct
{
  char **nn, *name, *seq;
  int n_nn, n_score, *score;
  size_t nchars;
};

struct cluster_struct
{
  fastaseq_t *fs;
  int n_fs, n_score, *idx, n_idx;
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

struct query_struct
{
  alignment aln;
  char *consensus;
  size_t *idx_c, *idx_m, *idx, trim; // idx_consensus (strictly same base present in all), idx_missing (consensus, but a gap/invalid in some sequence), idx_unique (must verify each query)
  int n_idx_c, n_idx_m, n_idx, dist;
  bool acgt;
};

int compare_fastaseq (const void *a, const void *b); /* \brief more neighours first, compare_score() if tie */
int compare_fastaseq_score (const void *a, const void *b);
fastaseq_t new_fastaseq (int n_score);
void del_fastaseq (fastaseq_t fs);
void update_fasta_seq (fastaseq_t to, char **seq, char **name, size_t nchars, int *score);

cluster_t new_cluster (char *seq, size_t nchars, int mindist, size_t trim, int n_score);
void del_cluster (cluster_t clus);
void generate_idx_from_cluster_list (cluster_t *clust, int n_clust, int min_freq);
void check_seq_against_cluster (cluster_t clust, char **seq, char **name, size_t nchars);
void add_seq_to_cluster (cluster_t clust, int idx, char **seq, char **name, size_t nchars, int *score);
int merge_clusters (cluster_t clust1, cluster_t clust2); 
int compact_cluster (cluster_t clust);

void save_cluster_to_xz_file (cluster_t *clust, int n_clust, const char* filename);
void save_cluster_to_gz_file (cluster_t *clust, int n_clust, const char* filename);
void save_neighbours_to_xz_file (cluster_t *clust, int n_clust, const char* filename);
void save_neighbours_to_gz_file (cluster_t *clust, int n_clust, const char* filename);

readfasta_t new_readfasta (const char *seqfilename);
int readfasta_next (readfasta_t rfas);
void del_readfasta (readfasta_t rfas);
int accumulate_reference_sequence (char **ref, char *s, size_t nsites);
int replace_Ns_from_reference (char *ref, size_t nsites);
void quick_pairwise_score_acgt_and_valid (char *s1, char *s2, size_t nsites, int maxdist, int *score, size_t *idx); // the others are used internally, tis the only one used upstream
int quick_count_sequence_non_N (char *s, size_t nsites);

/* for uvaia_ball */
void seq_ball_against_query_structure (char **seq, int *min_dist, int ball_radius, query_t qu);
query_t new_query_structure_from_fasta (char *filename, int trim, int dist, int acgt);
void del_query_structure (query_t qu);
void create_query_indices (query_t qu);
void reorder_query_structure (query_t qu);
void exclude_redundant_query_sequences (query_t qu, int keep_more_resolved);
#endif
