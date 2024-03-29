/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "utils.h"

void describe_scores (char *query_name, double *score, char_vector refnames, int nbest, int nmax, int **idx, int *n_idx);
bool sequence_n_below_threshold (char *seq, int seq_length, double threshold);

void
uvaia_keep_only_valid_sequences (alignment aln, double ambiguity, bool check_aligned)
{
  double result[3];
  int i, *valid = NULL, n_valid=0, maxlength = 0;
  valid = (int*) biomcmc_malloc (aln->ntax * sizeof (int));

  for (i=0; i < aln->character->nstrings; i++) {
    if (aln->character->nchars[i] < 5) {
      fprintf (stderr, "Sequence %s is too short ( = %lu sites), limit is hardcoded at 5bps.\n", aln->taxlabel->string[i], aln->character->nchars[i]);
      continue;
    }
    upper_kseq (aln->character->string[i], aln->character->nchars[i]); // must be _before_ count_sequence
    biomcmc_count_sequence_acgt (aln->character->string[i], aln->character->nchars[i], result); 
    if (result[2] > ambiguity) { 
      fprintf (stderr, "Sequence %s has proportion of N etc. (=%lf) above threshold of %lf\n", aln->taxlabel->string[i], result[2], ambiguity);
      continue;
    }
    if (result[0] < 1. - 1.1 * ambiguity) { 
      fprintf (stderr, "Sequence %s has proportion of ACGT (=%lf) below threshold of %lf\n", aln->taxlabel->string[i], result[0], 1. - 1.1 * ambiguity); 
      continue;
    }
    valid[n_valid++] = i;
    if (!maxlength) maxlength = aln->character->nchars[i];
    else if ((maxlength != -1) && ((size_t) maxlength != aln->character->nchars[i])) maxlength = -1;
  } // for each string
  
  /* check if remaining sequences have same length (if aln->is_aligned is true it implies patterns are used, which we don't here) */
  if (check_aligned && maxlength == -1) {
    biomcmc_warning ("Reference sequences in file %s are not aligned.\n", aln->filename);
    del_alignment (aln);
    biomcmc_error ("You can use uvaialign (or mafft, or minimap2) to align them against the same reference.");
  }

  /* alignment will not be fully functional since we dont update taxshort and taxlabel_hash  */
  char_vector_reduce_to_valid_strings (aln->character, valid, n_valid);
  char_vector_reduce_to_valid_strings (aln->taxlabel, valid, n_valid);
  aln->ntax = n_valid;
}

double
query_genome_against_char_vectors (char *name, char *s, unsigned l, char_vector cv_seq, char_vector cv_name, int nbest, int nmax, int **idx, int *n_idx, size_t trim)
{
  int i;
  int64_t time0[2];
  double *score;

  biomcmc_get_time (time0);
  score = (double*) biomcmc_malloc ((5 * cv_seq->nstrings) * sizeof (double)); // 3 scores but we need 5 for openMP

  if (cv_seq->nchars[0] != (size_t) l) { // all refs have same size
    biomcmc_warning ("this program assumes aligned sequences, and sequence %s has length %u while reference sequences have length %lu", name, l, cv_seq->nchars[0]);
    if (score) free (score);
    return biomcmc_update_elapsed_time (time0); // returns time in seconds
  }

#ifdef _OPENMP
#pragma omp parallel for shared(score, time0, cv_seq, cv_name) schedule(dynamic)
#endif
  for (i = 0; i < cv_seq->nstrings; i++) { // openmp doesnt like complex for() loops
    biomcmc_pairwise_score_matches (cv_seq->string[i] + trim, s + trim, l - 2 * trim, score + (5 * i));
//    score[5 * i] = result[0]; // ACGT match
    score[5 * i + 1] /= score[5 * i + 4]; // match between non-N (i.e. include R W etc.) divided by compared sites
//    score[5 * i + 2] = result[2]; // weighted partial matches (e.g. T has 50% match  with W (T+A)) 
//    if (score[5*i+2] < 1.) { for (int j=0; j < 5; j++) printf ("%lf ", score[5*i + j]); printf ("DBG:: %u \n", l); }
  }
  describe_scores (name, score, cv_name, nbest, nmax, idx, n_idx); // updates list idx[] of references to save
  if (score)  free (score);
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
  for (i=0; i < refnames->nstrings; i++) efd->d[i].freq = score[5 * i];
  sort_empfreq_double_decreasing (efd); 
  best_score = efd->d[0].freq;
  for (i = 0; i < nbest; i++) idbest[i] = efd->d[i].idx;
  for (; (i < nmax) && ((best_score  - efd->d[i].freq) < 1.e-15); i++) idbest[i] = efd->d[i].idx;
  n_idbest = i;

  /* 2. find references using score[1]/n_valid (non-N matches over non-N) excluding lower score[0] sequences */
  if (refnames->nstrings > 16) {
    empfreq_double efd2 = new_empfreq_double (refnames->nstrings/8);
    for (i=0; i < efd2->n; i++) {
      efd2->d[i].idx = efd->d[i].idx; // idx of seqs with best score
      efd2->d[i].freq = score[5 * efd->d[i].idx + 1]; // score[1] 
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
      efd2->d[i].freq = score[5 * efd->d[i].idx + 2]; // score[2] 
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

  /* 4. idbest[] will have some duplicate ids, we remove those and store in order of frequency (over scores) */
  empfreq bid = new_empfreq_from_int (idbest, n_idbest);
  /* 4.1 efd2 will reorder bid by ACGT_match score */ 
  empfreq_double efd2 = new_empfreq_double (bid->n);
  for (i=0; i < efd2->n; i++) {
    efd2->d[i].idx = bid->i[i].idx; // idx of seq to print  
    efd2->d[i].freq = score[5 * bid->i[i].idx]; // score[0] is ACGT_match 
  }
  sort_empfreq_double_decreasing (efd2);
  /* 4.2 order of efd2 is used to print */
  for (i = 0; i < bid->n; i ++) 
    printf ("%48s, %48s, %13.0lf, %13.1lf, %13.8lf, %13.3lf\n", query_name, refnames->string[efd2->d[i].idx],  score[5 * efd2->d[i].idx + 4],
            score[5 * efd2->d[i].idx], score[5 * efd2->d[i].idx + 1], score[5 * efd2->d[i].idx + 2]);
  fflush(stdout);

  (*idx) = (int*) biomcmc_realloc ((int*) (*idx), ((*n_idx) + bid->n) * sizeof (int));
  for (i = 0; i < bid->n; i++) (*idx)[(*n_idx)++] = bid->i[i].idx;
  del_empfreq_double (efd);
  del_empfreq_double (efd2);
  del_empfreq (bid);
  return;
}

void
print_score_header (void)
{    
  printf ("%48s, %48s, %13s, %13s, %13s, %13s\n","query sequence", "reference sequence", "valid_sites", "ACGT_matches", "prop_char_matches", "partial_matches");
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
  save_xzfasta_from_char_vector (filename, name, seq);

  del_empfreq (best_ids);
}

// below are unused/testing functions

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

char *
return_query_aligned (int pattern_length, char* text, int text_length, edit_cigar_t* edit_cigar, mm_allocator_t* mm_allocator) 
{
  // Parameters
  char* const operations = edit_cigar->operations;
  // Allocate alignment buffers
  const int max_buffer_length = text_length + pattern_length + 1;
  char* text_alg = mm_allocator_calloc (mm_allocator, max_buffer_length, char, true);
  // Compute alignment buffers
  int i, alg_pos = 0, pattern_pos = 0, text_pos = 0;
  for (i = edit_cigar->begin_offset; i < edit_cigar->end_offset; ++i) switch (operations[i]) {
    case 'M':
    case 'X':
      text_alg[alg_pos++] = text[text_pos++];
      pattern_pos++;
      break;
    case 'I':
      text_pos++;
      break;
    case 'D':
      text_alg[alg_pos++] = '-';
      pattern_pos++;
      break;
    default:
      break;
  }
  text_alg[alg_pos] = '\0';
  return text_alg;
}

static uint8_t is_acgt[256] = {0xff};
static uint8_t is_indel[256] = {0xff};

void
initialise_acgt (void)
{
  for (int i = 0; i < 256; i++) is_acgt[i] = is_indel[i] = 0;
  is_acgt['A'] = is_acgt['C'] = is_acgt['G'] = is_acgt['T'] = is_acgt['a'] = is_acgt['c'] = is_acgt['g'] = is_acgt['t'] = 1;
  is_indel['N'] = is_indel['n'] = is_indel['X'] = is_indel['x'] = is_indel['-'] = is_indel['?'] = is_indel['O'] = is_indel['o']  = is_indel['.'] = 1;
}

int
is_site_acgt_distinct_pair (char s1, char s2) // relies on external call to initialise_acgt()
{
  if ((!is_acgt[(int)s1]) || (!is_acgt[(int)s2])) return 0;
  return (s1 != s2); // assume toupper()
}

int
is_site_acgt_pair_valid (char s1, char s2) // relies on external call to initialise_acgt()
{
  return ((is_acgt[(int)s1]) && (is_acgt[(int)s2]));
}

int
is_site_pair_valid (char s1, char s2) // relies on external call to initialise_acgt()
{
  return ((!is_indel[(int)s1]) && (!is_indel[(int)s2]));
}

int
is_site_acgt (char s1) // relies on external call to initialise_acgt()
{
  return (int) is_acgt[(int)s1];
}

int
is_site_valid (char s1) // relies on external call to initialise_acgt()
{
  return (int) !is_indel[(int)s1];
}

void // not used (copied from first pilot "search_u.c" )
compress_kseq_to_acgt (char *s, unsigned *l)
{
  unsigned int i, new_l = 0;
  for (i = 0; i < *l; i++) if (is_acgt[ (int) s[i] ]) s[new_l++] = toupper(s[i]);
  *l = new_l;
  //  seq->seq.s = biomcmc_realloc ((char*)seq->seq.s, (new_l + 1) * sizeof (char)); // seq calloc size is seq.m not seq.l
  s[new_l] = '\0';
}

double 
calculate_score_from_cigar (edit_cigar_t* edit_cigar) 
{
  int i, matches = 0; 
  for (i = edit_cigar->begin_offset; i < edit_cigar->end_offset; ++i) if (edit_cigar->operations[i] == 'M') matches++; 
  return (double)(matches)/(double)(edit_cigar->end_offset - edit_cigar->begin_offset);
}
