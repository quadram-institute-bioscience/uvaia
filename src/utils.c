/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "utils.h"

void upper_kseq (char *s, unsigned l);
void describe_scores (char *query_name, double *score, char_vector refnames, int nbest, int nmax, int **idx, int *n_idx);
bool sequence_n_below_threshold (char *seq, int seq_length, double threshold);


void
add_reference_genome_to_char_vectors (char *name, char *s, unsigned l, char_vector cv_seq, char_vector cv_name)
{
  double result[3];
  if (cv_seq->next_avail && (cv_seq->nchars[cv_seq->next_avail-1] != (size_t) l)) {
    biomcmc_warning ("This program assumes aligned reference sequences, and sequence %s has length %u while sequence %s has length %lu\n",
                     name, l, cv_name->string[cv_seq->next_avail-1], cv_seq->nchars[cv_seq->next_avail-1]);
    biomcmc_error ("You can use uvaia_align (or mafft, or minimap2) to align them against the same reference");
  }
  upper_kseq (s, l); // must be _before_ count_sequence
  biomcmc_count_sequence_acgt (s, l, result); 
  if (result[0] < 0.5) { fprintf (stderr, "Reference %s has proportion of ACGT (=%lf) below 50%% threshold\n", name, result[0]); return; }
  if (result[2] > 0.5) { fprintf (stderr, "Reference %s has proportion of N etc. (=%lf) above 50%% threshold\n", name, result[2]); return; }
  char_vector_add_string (cv_seq,  s);
  char_vector_add_string (cv_name, name);
}

double
query_genome_against_char_vectors (char *name, char *s, unsigned l, char_vector cv_seq, char_vector cv_name, int nbest, int nmax, int **idx, int *n_idx)
{
  int i;
  int64_t time0[2];
  double *score, result[3];

  biomcmc_get_time (time0);
  upper_kseq (s, l); // must be _before_ count_sequence
  biomcmc_count_sequence_acgt (s, l, result); 
  if (result[0] < 0.5) { fprintf (stderr, "Query %s has proportion of ACGT (=%9lf) below 50%% threshold\n", name, result[0]); return 0.; }
  if (result[2] > 0.5) { fprintf (stderr, "Query %s has proportion of N etc. (=%9lf) above 50%% threshold\n", name, result[2]); return 0.; }
  score = (double*) biomcmc_malloc (3 * cv_seq->nstrings * sizeof (double)); // 3 scores

  if (cv_seq->nchars[0] != (size_t) l) { // all refs have same size
    biomcmc_warning ("this program assumes aligned sequences, and sequence %s has length %u while reference sequence %s has length %lu",
                     name, l, cv_name->string[0], cv_seq->nchars[0]);
    biomcmc_error ("Soon we'll be able to align ourselves ;)");
  }

#ifdef _OPENMP
#pragma omp parallel for shared(score, time0, cv_seq, cv_name) schedule(dynamic)
#endif
  for (i = 0; i < cv_seq->nstrings; i++) { // openmp doesnt like complex for() loops
    score[3 * i] = biomcmc_pairwise_score_matches (cv_seq->string[i], s, l, result); // ACGT match
    score[3 * i + 1] = result[0]; // match between non-N (i.e. include R W etc.)
    score[3 * i + 2] = result[1]; // unweighted partial matches (e.g. T matches partially with W (T+A) ) 
  }
  describe_scores (name, score, cv_name, nbest, nmax, idx, n_idx); // updates list idx[] of references to save
  if (score) free (score);
  return biomcmc_update_elapsed_time (time0); // returns time in seconds
}

void
upper_kseq (char *s, unsigned l)
{
  for (unsigned i = 0; i < l; i++) s[i] = toupper(s[i]);
}

void 
describe_scores (char *query_name, double *score, char_vector refnames, int nbest, int nmax, int **idx, int *n_idx)
{
  int i, *idbest, n_idbest;
  double best_score;

  if (nbest > refnames->nstrings) nbest = refnames->nstrings;
  if (nmax  > refnames->nstrings) nmax = refnames->nstrings;
  idbest = (int*) biomcmc_malloc (3 * nmax * sizeof (int)); 

  /* 1. find best references using score[0] (ACGT strict match) */
  empfreq_double efd = new_empfreq_double (refnames->nstrings);
  for (i=0; i < refnames->nstrings; i++) efd->d[i].freq = score[3 * i];
  sort_empfreq_double_decreasing (efd); 
  best_score = efd->d[0].freq;
  for (i = 0; i < nbest; i++) idbest[i] = efd->d[i].idx;
  for (; (i < nmax) && ((best_score  - efd->d[i].freq) < 1.e-15); i++) idbest[i] = efd->d[i].idx;
  n_idbest = i;

  /* 2. find references using score[1] (non-N matches) excluding lower score[0] sequences */
  if (refnames->nstrings > 16) {
    empfreq_double efd2 = new_empfreq_double (refnames->nstrings/8);
    for (i=0; i < efd2->n; i++) {
      efd2->d[i].idx = efd->d[i].idx; // idx of seqs with best score
      efd2->d[i].freq = score[3 * efd->d[i].idx + 1]; // score[1] 
    }
    sort_empfreq_double_decreasing (efd2);
    if (nbest > efd2->n) nbest = efd2->n;
    if (nmax  > efd2->n) nmax  = efd2->n;
    best_score = efd2->d[0].freq;
    for (i = 0; i < nbest; i++) idbest[i + n_idbest] = efd2->d[i].idx;
    for (; (i < nmax) && ((best_score  - efd2->d[i].freq) < 1.e-15); i++) idbest[i + n_idbest] = efd2->d[i].idx;
    n_idbest += i;
    del_empfreq_double (efd2);
  }

  /* 3. find references using score[2] (partial non-N matches) excluding lower score[0] sequences */
  if (refnames->nstrings > 32) {
    empfreq_double efd2 = new_empfreq_double (refnames->nstrings/16);
    for (i=0; i < efd2->n; i++) {
      efd2->d[i].idx = efd->d[i].idx; // idx of seqs with best score
      efd2->d[i].freq = score[3 * efd->d[i].idx + 2]; // score[2] 
    }
    sort_empfreq_double_decreasing (efd2);
    if (nbest > efd2->n) nbest = efd2->n;
    if (nmax  > efd2->n) nmax  = efd2->n;
    best_score = efd2->d[0].freq;
    for (i = 0; i < nbest; i++) idbest[i + n_idbest] = efd2->d[i].idx;
    for (; (i < nmax) && ((best_score  - efd2->d[i].freq) < 1.e-15); i++) idbest[i + n_idbest] = efd2->d[i].idx;
    n_idbest += i;
    del_empfreq_double (efd2);
  }

  /* 4. idbest[] will have some duplicate ids, we remove those and print in order of frequency (over scores) */
  empfreq bid = new_empfreq_from_int (idbest, n_idbest);
  for (i = 0; i < bid->n; i ++)
    printf ("%48s, %48s, %15lf %15lf %15lf\n", query_name, refnames->string[bid->i[i].idx], 
            score[3 * bid->i[i].idx], score[3 * bid->i[i].idx + 1], score[3 * bid->i[i].idx + 2]);
  fflush(stdout);

  (*idx) = (int*) biomcmc_realloc ((int*) (*idx), ((*n_idx) + bid->n) * sizeof (int));
  for (i = 0; i < bid->n; i++) (*idx)[(*n_idx)++] = bid->i[i].idx;
  del_empfreq_double (efd);
  del_empfreq (bid);
  return;
}

void
print_score_header (void)
{    
  printf ("%48s, %48s, %15s %15s %15s\n","query sequence", "reference sequence", "ACGT_matches", "char_matches", "partial_matches");
}

void 
describe_scores_old (char *query_name, double *score, char_vector refnames, int nbest, int nmax, int **idx, int *n_idx)
{
  int i; 
  empfreq_double srt = new_empfreq_double_sort_decreasing (score, refnames->nstrings);
  double best_score = srt->d[0].freq;
  if (nbest > refnames->nstrings) nbest = refnames->nstrings;
  if (nmax > refnames->nstrings)  nmax = refnames->nstrings;

  for (i = 0; i < nbest; i++) 
    printf ("%s, %s, %lf\n", query_name, refnames->string[srt->d[i].idx], srt->d[i].freq);
  for (; (i < nmax) && ((best_score  - srt->d[i].freq) < 1.e-15); i++) 
    printf ("%s, %s, %lf\n", query_name, refnames->string[srt->d[i].idx], srt->d[i].freq);
  fflush(stdout);
  nbest = i; // in case we had more than nbest with optimal score 

  (*idx) = (int*) biomcmc_realloc ((int*) (*idx), ((*n_idx) + nbest) * sizeof (int));
  for (i = 0; i < nbest; i++) (*idx)[(*n_idx)++] = srt->d[i].idx;
  del_empfreq_double (srt);
  return;
}

void
save_sequences (const char *filename, int *idx, int n_idx, char_vector seq, char_vector name)
{
  int i;
  empfreq best_ids = new_empfreq_from_int (idx, n_idx); // remove duplicates (ref sequences chosen by several queries)
  n_idx = best_ids->n;
  for (i = 0; i < n_idx; i ++) idx[i] = best_ids->i[i].idx; // most frequent first
  qsort (idx, n_idx, sizeof (int), compare_int_increasing); // char_vector_reduce() assumes ordered idx of valid
  char_vector_reduce_to_valid_strings (seq, idx, n_idx);
  char_vector_reduce_to_valid_strings (name, idx, n_idx);
  save_gzfasta_from_char_vector (filename, name, seq);

  del_empfreq (best_ids);
}

bool
sequence_n_below_threshold (char *seq, int seq_length, double threshold)
{
  int count = 0, i;
  char this;
  if (threshold < 1e-9) threshold = 1e-9;
  if (threshold > 1.)   threshold = 1.;
  int max_Ns = (int)(threshold * (double)(seq_length));

  for (i = 0; (i < seq_length) && (count < max_Ns); i++) {
    this = toupper(seq[i]);
    if ((this =='N') || (this == '-') || (this == '?') || (this == '.')) count++;
  }
  if (count < max_Ns) return true;
  return false;
}

