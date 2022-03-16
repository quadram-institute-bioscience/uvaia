/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * The approach used here for pointers is dangerous since the functions modify the pointer status of strings (the proper
 * usage involves reference_counters and delete()); i.e. one function allocs memory while another frees it as a
 * secondary effect. Strings are pointed to NULL to make sure they cannot modify the string anymore. 
 */

#include "fastaseq.h"


int compare_fastaseq (const void *a, const void *b);
int compare_fastaseq_score (const void *a, const void *b);
void quick_pairwise_score_reference (char *s1, char *s2, int nsites, int *score, int n_score, int *counter);
void quick_pairwise_score_truncated (char *s1, char *s2, int nsites, int maxdist, int *score);
void quick_pairwise_score_truncated_idx (char *s1, char *s2, size_t nsites, int maxdist, int *score, int *idx);
void quick_pairwise_score_truncated_idx_indelcheck (char *s1, char *s2, size_t nsites, int maxdist, int *score, size_t *idx);
void quick_pairwise_score_acgt (char *s1, char *s2, size_t nsites, int maxdist, int *score, size_t *idx);
int quick_count_sequence_non_N (char *s, size_t nsites);
int left_is_resolved_right (char *s1, char *s2, size_t nsites, size_t *idx);
int left_is_resolved_right_acgt (char *s1, char *s2, size_t nsites, size_t *idx);

int
compare_fastaseq (const void *a, const void *b)
{  
  int res_i = (*(fastaseq_t*)b)->n_nn - (*(fastaseq_t*)a)->n_nn; // decreasing since sort is ascending order (-1 to +1)
  if (res_i) return res_i;
  return compare_fastaseq_score (a,b);
}

int
compare_fastaseq_score (const void *a, const void *b)
{  
  int i, res;
  for (i = 0; i < (*(fastaseq_t*)b)->n_score + 2; i++) { // decreasing distance from reference, distance to first SNPs, non-Ns
    res = (*(fastaseq_t*)b)->score[i] - (*(fastaseq_t*)a)->score[i];
    if (res) return res;
  }
  return -1;
}

fastaseq_t
new_fastaseq (int n_score)
{
  int i;
  fastaseq_t fs = (fastaseq_t) biomcmc_malloc (sizeof (struct fastaseq_struct));
  fs->n_nn = 0;
  if (n_score < 0) n_score = 0;
  fs->n_score = n_score; // position of first SNPs 
  fs->score = (int*) biomcmc_malloc ((size_t)(n_score + 2) * sizeof (int)); // +1 for total distance (= score[0] )
  // score = [dist_to_ref, loc_of_SNP_1, loc_of_SNP_2,..., loc_of_SNP_n, number of ACGT+]
  for (i=0; i < n_score + 2; i++) fs->score[i] = 0;
  fs->nchars = 0;
  fs->name = fs->seq = NULL;
  fs->nn = NULL;
  return fs;
}

void
del_fastaseq (fastaseq_t fs)
{
  if (!fs) return;
  int i;
  if (fs->nn) {
    for (i = fs->n_nn-1; i >= 0; i--) if (fs->nn[i]) free (fs->nn[i]);
    free (fs->nn);
  }
  if (fs->score) free (fs->score);
  if (fs->name)  free (fs->name);
  if (fs->seq)   free (fs->seq);
  free (fs);
  return;
}

void
update_fasta_seq (fastaseq_t to, char **seq, char **name, size_t nchars, int *score)
{
  // 1. replace sequence info 
  if (to->seq) free (to->seq);
  to->seq = *seq;
  *seq = NULL;
  // 2. update list of neighbours
  if (to->name) {
    to->nn = (char**) biomcmc_realloc ((char**) to->nn, (to->n_nn + 1) * sizeof (char*));
    to->nn[to->n_nn++] = to->name;
  }
  // 3. replace name
  to->name = *name;
  *name = NULL;
  // 4. update score and seq length info
  for (int i = 0; i < to->n_score + 2; i++) to->score[i] = score[i];
  to->nchars = nchars; // we assume checks were done outside this function; we accept new sequence length
}

cluster_t
new_cluster (char *seq, size_t nchars, int mindist, size_t trim, int n_score)
{
  cluster_t clust = (cluster_t) biomcmc_malloc (sizeof (struct cluster_struct));
  clust->fs = NULL;
  clust->n_fs = 0;
  clust->mindist = mindist;
  clust->n_score = n_score;
  clust->trim = trim;
  clust->reference = (char*) biomcmc_malloc ((nchars+1) * sizeof (char));
  memcpy (clust->reference, seq, nchars);
  clust->reference[nchars] = '\0';
  clust->nchars = nchars;
  clust->n_idx = (int) nchars;
  clust->idx = (int*) biomcmc_malloc (nchars * sizeof (int)); 
  for (size_t i = 0; i < nchars; i++) clust->idx[i] = 0; // initially holds counts if SNPs wrt ref per site; later will have idx of those

  return clust;
}

void
del_cluster (cluster_t clust)
{
  if (!clust) return;
  int i;
  for (i = clust->n_fs-1; i >= 0; i--) del_fastaseq (clust->fs[i]);
  if (clust->reference) free (clust->reference);
  if (clust->idx) free (clust->idx);
  free (clust);
  return;
}

void
generate_idx_from_cluster_list (cluster_t *clust, int n_clust, int min_freq)
{
  int c, i, n_i = 0;
  if (min_freq < 0) min_freq = 0;
  for (c = 1; c < n_clust; c++) for (i = 0; i < clust[0]->n_idx; i++) clust[0]->idx[i] += clust[c]->idx[i];
  for (i = 0; i < clust[0]->n_idx; i++) if (clust[0]->idx[i] > min_freq) clust[0]->idx[n_i++] = i;
  for (c = 0; c < n_clust; c++) clust[c]->n_idx = n_i;
  for (c = 0; c < n_clust; c++) clust[c]->idx = (int*) biomcmc_realloc ((int*) clust[c]->idx, (size_t) (n_i) * sizeof (int)); 
  for (c = 1; c < n_clust; c++) for (i = 0; i < clust[0]->n_idx; i++) clust[c]->idx[i] = clust[0]->idx[i];
  return;
}

void
check_seq_against_cluster (cluster_t clust, char **seq, char **name, size_t nchars)
{
  int i, minloc = 0;
  size_t scorelength = (size_t)(clust->n_score) + 2;
  int *score; // score = {dist to reference, locations of first SNPs w.r.t. ref, number of ACGT (non-N actually)}

  if (nchars != clust->nchars) 
    biomcmc_error ("%s cannot work with unaligned sequences; sequence %s has %lu sites while reference has %d.", PACKAGE_STRING, *name, nchars, clust->nchars);

  score = (int*) biomcmc_malloc (scorelength * sizeof (int));
  upper_kseq (*seq, nchars); // must be _before_ count_sequence
  score[scorelength-1] = quick_count_sequence_non_N (*seq, nchars); // last element of score vector is not touched by quick_pairwise (unless you provide a larger n_score)
  quick_pairwise_score_reference (*seq + clust->trim, clust->reference + clust->trim, (int) nchars - 2 * clust->trim, score, clust->n_score, clust->idx);
  for (i = 0; i < clust->n_fs; i++) if (abs(score[0] - clust->fs[i]->score[0]) <= clust->mindist) { // only if within ball ring from reference
    if (clust->n_score) minloc = BIOMCMC_MIN (score[1], clust->fs[i]->score[1]) - 1; // location of first disagreement with reference
    else minloc = 0;
    if (minloc < 0) minloc = 0;
    quick_pairwise_score_truncated (*seq + clust->trim + minloc, clust->fs[i]->seq + clust->trim + minloc, (int) nchars - 2 * clust->trim, clust->mindist+1, score);
    if (score[0] <= clust->mindist) { // text matches (hamming distance) (BTW this distance score[0] is not used downstream) 
      if (score[0]) score[scorelength] = 0.;  // seq has SNP wrt current medoid; pretend it's all N so just name (but not sequence) is stored in neighbour list
      add_seq_to_cluster (clust, i, seq, name, nchars, score);
      if (score) free (score);
      return;
    }
  }
  // not similar to any existing cluster
  add_seq_to_cluster (clust, i, seq, name, nchars, score);
  if (score) free (score);
  return;
}

void
add_seq_to_cluster (cluster_t clust, int idx, char **seq, char **name, size_t nchars, int *score)
{
  if (idx >= clust->n_fs) { // new cluster, no similar sequences in cluster set
    idx = clust->n_fs;
    clust->fs = (fastaseq_t*) biomcmc_realloc ((fastaseq_t*) clust->fs, (++clust->n_fs) * sizeof (fastaseq_t));
    clust->fs[idx] = new_fastaseq (clust->n_score);
    update_fasta_seq (clust->fs[idx], seq, name, nchars, score);
    return;
  }
  if (score[clust->n_score+1] > clust->fs[idx]->score[clust->n_score+1]) { // existing cluster, this sequence is new medoid since furthest from reference
    update_fasta_seq (clust->fs[idx], seq, name, nchars, score);
    return;
  }
  // existing cluster and current sequence not new medoid
  clust->fs[idx]->nn = (char**) biomcmc_realloc ((char**) clust->fs[idx]->nn, (clust->fs[idx]->n_nn + 1) * sizeof (char*));
  clust->fs[idx]->nn[clust->fs[idx]->n_nn++] = *name;
  *name = NULL;
  if (*seq) free (*seq);
  *seq = NULL;
  return;
}

int
merge_clusters (cluster_t clust1, cluster_t clust2)
{
  int i, j, k, *score, count = 0, first = 0, last = 0, c2s, c1_n_fs = clust1->n_fs, *dst, maxdst, *idx2;
  fastaseq_t f1, f2;

  qsort (clust1->fs, clust1->n_fs, sizeof (fastaseq_t), compare_fastaseq_score); // decreasing score (i.e. distance from reference)
  qsort (clust2->fs, clust2->n_fs, sizeof (fastaseq_t), compare_fastaseq_score); // cannot count on order since adding clust2 may change medoid dist to ref
  score = (int*) biomcmc_malloc ((size_t)(clust1->n_score + 2) * sizeof (int));

  /* preprocess list of dists indices for clust1 */
  maxdst = clust1->fs[0]->score[0]; // maximum distance to reference
  dst = (int*) biomcmc_malloc ((maxdst + 1) * sizeof (int)); // dst[i] = first index of distance i
  for (i = 0; i <= clust1->fs[0]->score[0]; i++) dst[i] = -1;
  for (i = 0; i < c1_n_fs; i++) if (dst[ clust1->fs[i]->score[0] ] < 0) dst[ clust1->fs[i]->score[0] ]  = i;
 
  /* preprocess list of dists indices for clust2 */
  idx2 = (int*) biomcmc_malloc (2 * (clust2->fs[0]->score[0] + 1) * sizeof (int)); // idx in cluster 1 of lower and upper dist bounds
  for (i = 0; i < 2 * (clust2->fs[0]->score[0] + 1); i++) idx2[i] = -1;
  for (j = 0; j < clust2->n_fs; j++) if (idx2[2 * clust2->fs[j]->score[0]] < 0) {
    c2s = clust2->fs[j]->score[0];
    for(i = c2s + clust1->mindist; (i >= 0) && (i <= maxdst) && (dst[i] < 0); i++); // as i increases dst[i] decreases
    if ((i >= 0) && (i <= maxdst)) first = dst[i]; // we have to care for both borders since clust1 and clust2 may not overlap... 
    else first = 0;                                 // e.g. clust1 dists = {0 ... 10} but clust2 dist = 15
    for(i = c2s - clust1->mindist - 1; (i >= 0) && (i <= maxdst) && (dst[i] < 0); i--); // as i decreases dst[i] increases; minus one since we want last idx before this distance
    if ((i >= 0) && (i <= maxdst)) last = dst[i];
    else last = c1_n_fs;
    idx2[2 * c2s] = first; 
    idx2[(2 * c2s) + 1] = last; 
    //printf ("DEBUG:: j=%d dist=%d first=%d last=%d dst[0]=%d dst[%d]=%d\n", j, c2s, first, last, dst[0], maxdst, dst[maxdst]);
  }
  //for (j=0; j <= clust2->fs[0]->score[0]; j++) printf ("idx2[%d] = [%3d %3d]\n", j, idx2[2*j], idx2[(2*j)+1]);
  
  for (j = 0; j < clust2->n_fs; j++) { 
    c2s = clust2->fs[j]->score[0];
    first = idx2[2 * c2s]; last = idx2[(2 * c2s) + 1];
    // c1_n_fs does not change when new clust2 elements are added to clust1 (no need to compare clust2 with itself)
    for (i = first; i < last; i++) if (abs(c2s - clust1->fs[i]->score[0]) <= clust1->mindist) { // only if within ball ring from reference 
      f1 = clust1->fs[i]; f2 = clust2->fs[j]; // nicknames to avoid silly mistakes
      quick_pairwise_score_truncated_idx (f1->seq + clust1->trim, f2->seq + clust1->trim, clust1->n_idx, clust1->mindist+1, score, clust1->idx);
      if (score[0] <= clust1->mindist) { // medoid will be decided by score[n_score+1] which was already calculated with or without reference
        add_seq_to_cluster (clust1, i, &(f2->seq), &(f2->name), f2->nchars, f2->score);
        if (f2->nn) {
          f1->nn = (char**) biomcmc_realloc ((char**) f1->nn, (f1->n_nn + f2->n_nn) * sizeof (char*));
          for (k = 0; k < f2->n_nn; k++) f1->nn[f1->n_nn + k] = f2->nn[k];
          f1->n_nn += f2->n_nn;
          free (f2->nn); f2->nn = NULL;  // avoid nn vector elements being freed by del_fastaseq() 
        }
        del_fastaseq (clust2->fs[j]); clust2->fs[j] = NULL;
        count++;
        break; // breaks clust1 loop 
      } // within distance
    } // for i in clust1
    if (i == last) { // clust2 not similar to any existing clust1; will increase clust1 size 
      clust1->fs = (fastaseq_t*) biomcmc_realloc ((fastaseq_t*) clust1->fs, (clust1->n_fs+1) * sizeof (fastaseq_t));
      clust1->fs[clust1->n_fs++] = clust2->fs[j];
      clust2->fs[j] = NULL;
    }
  } // for j in clust2  
  if (clust2->fs) free (clust2->fs);
  clust2->fs = NULL;
  clust2->n_fs = 0;
  if (score) free (score);
  if (dst) free (dst);
  if (idx2) free (idx2);
  return count;
}

int
compact_cluster (cluster_t clust) // not useful
{
  int i, j, k, *score, count = 0;
  fastaseq_t f1, f2;

  score = (int*) biomcmc_malloc ((size_t)(clust->n_score + 2) * sizeof (int));

  for (i = 0; i < clust->n_fs-1; i++) for (j = i + 1; j < clust->n_fs; j++) 
    if ((clust->fs[i] && clust->fs[j]) && (abs(clust->fs[i]->score[0] - clust->fs[j]->score[0]) <= clust->mindist)) { // only if within ball ring from reference 
    f1 = clust->fs[i]; f2 = clust->fs[j];
    quick_pairwise_score_truncated (f1->seq + clust->trim, f2->seq + clust->trim, (int) (clust->nchars - 2 * clust->trim), clust->mindist+1, score);
    if (score[0] <= clust->mindist) { 
      add_seq_to_cluster (clust, i, &(f2->seq), &(f2->name), f2->nchars, f2->score);
      if (f2->nn) {
        f1->nn = (char**) biomcmc_realloc ((char**) f1->nn, (f1->n_nn + f2->n_nn) * sizeof (char*));
        for (k = 0; k < f2->n_nn; k++) f1->nn[f1->n_nn + k] = f2->nn[k];
        f1->n_nn += f2->n_nn;
        free (f2->nn); f2->nn = NULL;  // avoid nn vector elements being freed by del_fastaseq() 
      }
      del_fastaseq (clust->fs[j]); clust->fs[j] = NULL;
      count++;
    } // within distance
  } // for i in clust1 and  for j in clust2  

  for (j = 0, i = 0; i < clust->n_fs; i++) if (clust->fs[i]) clust->fs[j++] = clust->fs[i]; // only non-null
  clust->n_fs = j;
  clust->fs = (fastaseq_t*) biomcmc_realloc ((fastaseq_t*) clust->fs, (clust->n_fs) * sizeof (fastaseq_t));
  if (score) free (score);
  return count;
}

void
save_cluster_to_xz_file (cluster_t *clust, int n_clust, const char* filename)
{
#ifndef HAVE_LZMA 
  fprintf (stderr, "LZMA library missing; reverting to gzip or uncompressed\n");
  save_cluster_to_gz_file (clust, n_clust, filename);
  return;
#else
  int i, c, errors = 0;
  size_t nchars;
  xz_file_t *xz = NULL;
  xz = biomcmc_xz_open (filename, "w", 4096);
  if (!xz) {
    fprintf (stderr, "Problem opening file %s for writing with XZ compression; reverting to gzip or uncompressed\n", filename);
    save_cluster_to_gz_file (clust, n_clust, filename);
    return;
  }
  for (c = 0; c < n_clust; c++) for (i = 0; i < clust[c]->n_fs; i++) {
    if (biomcmc_xz_write (xz, ">", 1) != 1) errors++; 
    nchars = strlen(clust[c]->fs[i]->name);
    if (biomcmc_xz_write (xz, clust[c]->fs[i]->name, nchars) != nchars) errors++;
    if (biomcmc_xz_write (xz, "\n", 1) != 1) errors++;
    nchars = clust[c]->fs[i]->nchars;
    if (biomcmc_xz_write (xz, clust[c]->fs[i]->seq, nchars) != nchars) errors++;
    if (biomcmc_xz_write (xz, "\n", 1) != 1) errors++;
  }
  biomcmc_xz_close (xz);
  if (errors) fprintf (stderr,"File %s may not be correctly compressed, %d error%s occurred.\n", filename, errors, ((errors > 1)? "s":""));
#endif
  return;
}

void
save_cluster_to_gz_file (cluster_t *clust, int n_clust, const char* filename)
{
  int i, c;
#ifdef HAVE_ZLIB
  gzFile stream;
  stream = biomcmc_gzopen (filename, "w");
  for (c = 0; c < n_clust; c++) for (i = 0; i < clust[c]->n_fs; i++) gzprintf (stream, ">%s\n%s\n", clust[c]->fs[i]->name, clust[c]->fs[i]->seq);
  gzclose (stream);
#else
  FILE *stream;
  stream = biomcmc_fopen (filename, "w");
  for (c = 0; c < n_clust; c++) for (i = 0; i < clust[c]->n_fs; i++) fprintf (stream, ">%s\n%s\n", clust[c]->fs[i]->name, clust[c]->fs[i]->seq);
  fclose (stream);
#endif
  return;
}

void
save_neighbours_to_xz_file (cluster_t *clust, int n_clust, const char* filename)
{
#ifndef HAVE_LZMA 
  fprintf (stderr, "LZMA library missing; reverting to gzip or uncompressed\n");
  save_neigbours_to_gz_file (clust, n_clust, filename);
  return;
#else
  int i, j, c, errors = 0;
  size_t nchars;
  xz_file_t *xz = NULL;
  xz = biomcmc_xz_open (filename, "w", 4096);
  if (!xz) {
    fprintf (stderr, "Problem opening neighbours file %s for writing with XZ compression; reverting to gzip or uncompressed\n", filename);
    save_neighbours_to_gz_file (clust, n_clust, filename);
    return;
  }
  for (c = 0; c < n_clust; c++) {
    for (i = 0; i < clust[c]->n_fs; i++) {
      nchars = strlen(clust[c]->fs[i]->name);
      if (biomcmc_xz_write (xz, clust[c]->fs[i]->name, nchars) != nchars) errors++; 
      for (j = 0; j < clust[c]->fs[i]->n_nn; j++) {
        if (biomcmc_xz_write (xz, ",", 1) != 1) errors++; 
        nchars = strlen(clust[c]->fs[i]->nn[j]);
        if (biomcmc_xz_write (xz, clust[c]->fs[i]->nn[j], nchars) != nchars) errors++;
      }
      if (biomcmc_xz_write (xz, "\n", 1) != 1) errors++; 
    }
  } // for cluster
  if (errors) fprintf (stderr,"File %s may not be correctly compressed, %d error%s occurred.\n", filename, errors, ((errors > 1)? "s":""));
  biomcmc_xz_close (xz);
#endif
  return;
}

void
save_neighbours_to_gz_file (cluster_t *clust, int n_clust, const char* filename)
{
  int i, j, c;
#ifdef HAVE_ZLIB
  gzFile stream;
  stream = biomcmc_gzopen (filename, "w");
  for (c = 0; c < n_clust; c++) {
    for (i = 0; i < clust[c]->n_fs; i++) {
      gzprintf (stream, "%s", clust[c]->fs[i]->name);
      for (j = 0; j < clust[c]->fs[i]->n_nn; j++) gzprintf (stream, ",%s", clust[c]->fs[i]->nn[j]);
      gzprintf (stream, "\n");
    }
  }
  gzclose (stream);
#else
  FILE *stream;
  stream = biomcmc_fopen (filename, "w");
  for (c = 0; c < n_clust; c++) {
    for (i = 0; i < clust[c]->n_fs; i++) {
      fprintf (stream, "%s", clust[c]->fs[i]->name);
      for (j = 0; j < clust[c]->fs[i]->n_nn; j++) fprintf (stream, ",%s", clust[c]->fs[i]->nn[j]);
      fprintf (stream, "\n");
    }
  }
  fclose (stream);
#endif
  return;
}

readfasta_t
new_readfasta (const char *seqfilename)
{
  readfasta_t rfas = (readfasta_t) biomcmc_malloc (sizeof (struct readfasta_struct));
  rfas->seqfile = biomcmc_open_compress (seqfilename, "r");
  rfas->next_name = rfas->name = rfas->seq = NULL;
  rfas->line_read = NULL;
  rfas->linelength = rfas->seqlength = 0;
  rfas->newseq = false;
  return rfas;
}

int
readfasta_next (readfasta_t rfas)
{
  size_t nchars = 0;
  char *line = NULL, *delim = NULL;
  if (!rfas->seqfile) return -1; // we've finished reading the file
  /* keep reading file */
  while (biomcmc_getline_compress (&(rfas->line_read), &(rfas->linelength), rfas->seqfile) != -1) {
    line = rfas->line_read; /* the variable *line_read should point always to the same value (no line++ or alike) */
    if (nonempty_fasta_line (line)) { 
      if ((delim = strchr (line, '>'))) { // new sequence; store its name and return previous sequence
        delim++; // skip symbol
        rfas->newseq = true; // next iteration will be new sequence (we're finished reading the one in buffer) 
        if (rfas->next_name) { // not the first sequence; thus we return name and seq
          if (rfas->name) free (rfas->name);
          rfas->name = rfas->next_name;
        }
        nchars = strlen (delim);
        rfas->next_name = (char*) biomcmc_malloc ((nchars + 1) * sizeof (char));
        memcpy (rfas->next_name, delim, nchars);
        rfas->next_name[nchars] = '\0';
        if (rfas->seqlength) return (int)(rfas->seqlength); // there is a seq in the buffer
      }
      else { // reading the sequence, which may be multi-line
        if (rfas->newseq) { // previous line was seq name, we are starting new one
          if (rfas->seq) free (rfas->seq);
          rfas->seq = NULL;
          rfas->seqlength = 0;
        }
        line = remove_space_from_string (line);
        line = uppercase_string (line);
        nchars = strlen (line);
        rfas->seq = (char*) biomcmc_realloc ((char*)rfas->seq, (rfas->seqlength + nchars + 1) * sizeof (char));
        memcpy (rfas->seq + rfas->seqlength, line, nchars); 
        rfas->seqlength += nchars;
        rfas->seq[rfas->seqlength] = '\0'; // notice that strlen is seqlength + nchars _plus_one_ (to fit the null char)
        rfas->newseq = false; // to make sure seq is not deleted 
      }
    } // nonempty()
  } // getline()
  // reached EOF; first prepare to return -1 on subsequent calls
  biomcmc_close_compress (rfas->seqfile);
  rfas->seqfile = NULL;
  if (rfas->line_read) free (rfas->line_read);
  rfas->line_read = NULL;
  // and then return whatever is in buffer (hopefully name and seq)
  if (rfas->next_name) { // not the first sequence; thus we return name and seq
    if (rfas->name) free (rfas->name);
    rfas->name = rfas->next_name;
    rfas->next_name = NULL; // otherwise will be pointing to the last name
  }
  return (int)(rfas->seqlength); // there should be a seq in the buffer
}

void
del_readfasta (readfasta_t rfas)
{
  if (!rfas) return;
  if (rfas->seqfile) biomcmc_close_compress (rfas->seqfile);
  if (rfas->line_read) free (rfas->line_read);
  if (rfas->seq) free (rfas->seq);
  if (rfas->next_name) free (rfas->next_name); // rfas->name is always just a pointer; was allocated with next_name
  if (rfas->name) free (rfas->name);
  free (rfas);
}

int
accumulate_reference_sequence (char **ref, char *s, size_t nsites)
{ // quick-and-dirty version of biomcmc would do
  size_t i=0;
  int count = 0; // number of non_ACGT

  if (!(*ref)) {
    *ref = biomcmc_malloc (nsites * sizeof (char));
    memcpy (*ref, s, nsites);
    for (i = 0; i < nsites; i++) {
      (*ref)[i] = toupper ((*ref)[i]);
      if (((*ref)[i] != 'A') && ((*ref)[i] != 'C') && ((*ref)[i] != 'G') && ((*ref)[i] != 'T')) { (*ref)[i] = 'N'; count++; }
    }
    return count;
  }

  for (i = 0; i < nsites; i++) {
    s[i] = toupper (s[i]);
    if ((*ref)[i] == 'N') {
      if ((s[i] == 'A') || (s[i] == 'C') || (s[i] == 'G') || (s[i] == 'T')) (*ref)[i] = s[i];
      else count++;
    }
  }
  return count;
}

int
replace_Ns_from_reference (char *ref, size_t nsites)
{
  int count = 0;
  for (size_t i = 0; i < nsites; i++) if (ref[i] == 'N') { ref[i] = 'A'; count++; }
  return count;
}

void
quick_pairwise_score_reference (char *s1, char *s2, int nsites, int *score, int n_score, int *counter)
{ // assumes upper(), and just count text matches; adds non-N to idx 
  int i;
  score[0] = 0;
  for (i = 1; i <= n_score; i++) score[i] = -1; 

  for (i=0; i < nsites; i++) {
    if ((s1[i] == 'N') || (s2[i] == 'N') || (s1[i] == '-') || (s2[i] == '-') || (s1[i] == '?') || (s2[i] == '?')) continue; // skip this column
    score[0]++;
    if (s1[i] == s2[i]) score[0]--;
    else counter[i]++; // SNP wrt reference (i.e. valid pair, which disagrees
    if ((n_score) && (score[0]) && (score[0] <= n_score) && (score[score[0]] < 0)) score[ score[0] ] = i; // n-th difference between s1 and s2
  }
  return;
}

void
quick_pairwise_score_truncated (char *s1, char *s2, int nsites, int maxdist, int *score)
{ // assumes upper(), and just count text matches; "score" actually means distance (the lower the better); skips indels
  int i;
  score[0] = 0;

  for (i=0; (i < nsites) && (score[0] < maxdist); i++) {
    if ((s1[i] == 'N') || (s2[i] == 'N') || (s1[i] == '-') || (s2[i] == '-') || (s1[i] == '?') || (s2[i] == '?')) continue; // skip this column
    score[0]++;
    if (s1[i] == s2[i]) score[0]--;
  }
  return;
}

void
quick_pairwise_score_truncated_idx (char *s1, char *s2, size_t nsites, int maxdist, int *score, int *idx)
{ // assumes upper(), and just count text matches;
  size_t i;
  *score = 0;
  for (i=0; (i < nsites) && (score[0] < maxdist); i++) if (s1[ idx[i] ] != s2[ idx[i] ]) (*score)++;
  return;
}

void
quick_pairwise_score_truncated_idx_indelcheck (char *s1, char *s2, size_t nsites, int maxdist, int *score, size_t *idx)
{ // assumes upper(), and just count text matches; unlike above, checks for indels
  size_t i, j;
  score[0] = 0;
  for (j=0; (j < nsites) && (score[0] < maxdist); j++) {
    i = idx[j];
    if ((s1[i] == 'N') || (s2[i] == 'N') || (s1[i] == '-') || (s2[i] == '-') || (s1[i] == '?') || (s2[i] == '?')) continue; // skip this column
    score[0]++;
    if (s1[i] == s2[i]) score[0]--;
  }
  return;
}

void
quick_pairwise_score_acgt (char *s1, char *s2, size_t nsites, int maxdist, int *score, size_t *idx)
{ // assumes upper(), and just count acgt matches; relies on call to initialise_acgt() from utils.c
  size_t j;
  score[0] = 0;
  for (j=0; (j < nsites) && (score[0] < maxdist); j++) score[0] += (int) is_site_acgt_distinct_pair (s1[idx[j]], s2[idx[j]]);
  return;
}

int
left_is_resolved_right (char *s1, char *s2, size_t nsites, size_t *idx) 
{ // minus one if left is resolved, plus one if right is resolved, zero if identical and larger than one if distinct
  size_t j;
  int score = 0;
  bool snp1, snp2;
  for (j=0; (j < nsites) && (score < 0xff); j++) {
    snp1 = !((s1[idx[j]] == 'N') || (s1[idx[j]] == '-') || (s1[idx[j]] == '?')); // snp1 = 1 if not a space; 0 o.w.
    snp2 = !((s2[idx[j]] == 'N') || (s2[idx[j]] == '-') || (s2[idx[j]] == '?'));
    if (snp1 == snp2) continue; // both indels or both SNPs (we assume same SNP, if called properly after calculating distance)
    if (snp1 > snp2) { 
      if (score > 0) return 0xff; // s2 was more resolved, i.e. more ACGTs, but here the ACGT is on s1 (so they're not equivalent)
      else score = -1;
    }
    if (snp1 < snp2) {
      if (score < 0) return 0xff; // s1 was more resolved, i.e. more ACGTs, but here the ACGT is on s2 (so they're not equivalent)
      else score = 1;
    }
  }
  return score;
}

int
left_is_resolved_right_acgt (char *s1, char *s2, size_t nsites, size_t *idx) 
{ // minus one if left is resolved, plus one if right is resolved, zero if identical and larger than one if distinct
  size_t j;
  int score = 0;
  bool snp1, snp2;
  for (j=0; (j < nsites) && (score < 0xff); j++) {
    snp1 = is_site_acgt(s1[idx[j]]);
    snp2 = is_site_acgt(s2[idx[j]]);
    if (snp1 == snp2) continue; // both indels or both SNPs (we assume same SNP, if called properly after calculating distance)
    if (snp1 > snp2) { 
      if (score > 0) return 0xff; // s2 was more resolved, i.e. more ACGTs, but here the ACGT is on s1 (so they're not equivalent)
      else score = -1;
    }
    if (snp1 < snp2) {
      if (score < 0) return 0xff; // s1 was more resolved, i.e. more ACGTs, but here the ACGT is on s2 (so they're not equivalent)
      else score = 1;
    }
  }
  return score;
}

int
quick_count_sequence_non_N (char *s, size_t nsites)
{
  int non_n = nsites;
  for (size_t i = 0; i < nsites; i++) if ((s[i] == 'N') || (s[i] == '-') || (s[i] == '?') || (s[i] == 'X') || (s[i] == 'O') || (s[i] == '.')) non_n--;
  return non_n;
}

int
quick_count_sequence_acgt (char *s, size_t nsites)
{
  int non_n = 0;
  for (size_t i = 0; i < nsites; i++) non_n += (int) is_site_acgt (s[i]);
  return non_n;
}

/* uvaia_ball */ 

void
seq_ball_against_query_structure (char **seq, int *min_dist, int ball_radius, query_t qu)
{
  int i, c_dist=0;
  if (qu->acgt) {
    quick_pairwise_score_acgt (*seq, qu->consensus, qu->n_idx_c, ball_radius, min_dist, qu->idx_c);
    if (*min_dist >= ball_radius) return;
    c_dist = *min_dist;
    for (i = 0; (i < qu->aln->ntax) && ((*min_dist + c_dist) >= ball_radius); i++) { 
      quick_pairwise_score_acgt (*seq, qu->aln->character->string[i], qu->n_idx, ball_radius - c_dist, min_dist, qu->idx);
    }
    *min_dist += c_dist;
    return;
  }
  else {
    quick_pairwise_score_truncated_idx_indelcheck (*seq, qu->consensus, qu->n_idx_c, ball_radius, min_dist, qu->idx_c);
    if (*min_dist >= ball_radius) return;
    c_dist = *min_dist;
    for (i = 0; (i < qu->aln->ntax) && ((*min_dist + c_dist) >= ball_radius); i++) { 
      quick_pairwise_score_truncated_idx_indelcheck (*seq, qu->aln->character->string[i], qu->n_idx, ball_radius - c_dist, min_dist, qu->idx);
    }
    *min_dist += c_dist;
    return;
  }
}

query_t
new_query_structure_from_fasta (char *filename, int trim, int dist, int acgt)
{
 query_t qu = (query_t) biomcmc_malloc (sizeof (struct query_struct));
  qu->consensus = NULL;
  qu->idx_c = qu->idx = NULL;
  qu->n_idx_c = qu->n_idx = 0;
  qu->acgt = (bool) acgt;
 
  qu->aln = read_fasta_alignment_from_file (filename, 0xf); // 0xf -> neither true or false, but gets only bare alignment info
  //int i,j;
  //for (i=0;i<30;i++) {for (j=0;j<30;j++) printf("%c",qu->aln->character->string[i][j]); printf ("\n"); }
  
  /* check if parameters (trim, min_dist) are compatible with sequence lenght */
  if (trim < 0) trim = 0; // lower bound is zero (default value); params.trim is an int but we need size_t (cq->trim)
  if (trim > qu->aln->nchar / 2.1) trim = qu->aln->nchar / 2.1; // if we trim more than 1/2 the genome there's nothing left
  qu->trim = (size_t) trim;
  if (dist < 0) dist = 0;
  if (dist > (qu->aln->nchar - 2 * trim)/10) dist = (int)((qu->aln->nchar - 2*trim)/10); 
  qu->dist = dist;

  return qu;
}

void
del_query_structure (query_t qu)
{
  if (!qu) return;
  if (qu->consensus) free (qu->consensus);
  if (qu->idx_c) free (qu->idx_c);
  if (qu->idx) free (qu->idx);
  del_alignment (qu->aln);
  free (qu);
}

void
create_query_indices (query_t qu)
{
  int i,j;
  char s2;
  initialise_acgt (); // from utils.c, has local vector 

  if (!qu->consensus) qu->consensus = (char*) biomcmc_malloc (qu->aln->nchar * sizeof (char)); // we may call this function more than once
  for (i = 0; i < qu->aln->nchar; i++) qu->consensus[i] = 'N';

  if (qu->acgt) for (i = (int) qu->trim; i < (qu->aln->nchar - (int) qu->trim); i++) for (j = 0; (j < qu->aln->ntax) && (qu->consensus[i] != '#'); j++) {
    s2 = qu->aln->character->string[j][i];
    if (!is_site_acgt (s2)) continue; // skip indels but also partially ambiguous, focusing on ACGT differences
    if (qu->consensus[i] == 'N') qu->consensus[i] = s2; // consensus is missing info, will receive ACGT from alignment
    else if (qu->consensus[i] != s2) qu->consensus[i] = '#'; // ACGT info from consensus conflicts with current sequence
  }
  else for (i = (int) qu->trim; i < (qu->aln->nchar - (int) qu->trim); i++) for (j = 0; (j < qu->aln->ntax) && (qu->consensus[i] != '#'); j++) {
    s2 = qu->aln->character->string[j][i];
    if ( (s2 == 'N') || (s2 == '-') || (s2 == '?')) continue; // skip this column (only difference from loop above, acgt)
    if (qu->consensus[i] == 'N') qu->consensus[i] = s2; // consensus is missing info, will receive wathever state is from alignment
    else if (qu->consensus[i] != s2) qu->consensus[i] = '#';
  }

  qu->idx_c = (size_t*) biomcmc_realloc ((size_t*) qu->idx_c, (qu->aln->nchar - (int) qu->trim) * sizeof (size_t)); 
  qu->idx   = (size_t*) biomcmc_realloc ((size_t*) qu->idx,   (qu->aln->nchar - (int) qu->trim) * sizeof (size_t)); 
  qu->n_idx_c = qu->n_idx = 0;
  
  for (i = (int) qu->trim; i < (qu->aln->nchar - (int) qu->trim); i++) if (qu->consensus[i] != 'N') {
    if (qu->consensus[i] == '#') qu->idx[qu->n_idx++] = i; // polymorphic site; we need to check every alignment sequence
    else                     qu->idx_c[qu->n_idx_c++] = i; // non-segregating site; enough to compare with consensus
  }
  qu->idx_c = (size_t*) biomcmc_realloc ((size_t*) qu->idx_c, qu->n_idx_c * sizeof (size_t)); 
  qu->idx   = (size_t*) biomcmc_realloc ((size_t*) qu->idx,   qu->n_idx   * sizeof (size_t)); 
}

void
reorder_query_structure (query_t qu)
{ // could use indices but this function is linear not quadratic on seqlength thus not critical
  int i, *non_n; // non_n can be number of non-N or number of ACGT; and later used to store indices of best-to-worse 
  empfreq e;

  non_n = (int*) biomcmc_malloc (qu->aln->ntax * sizeof (int));	
  if (qu->acgt) for (i = 0; i < qu->aln->ntax; i++) non_n[i] = quick_count_sequence_acgt  (qu->aln->character->string[i] + qu->trim, qu->aln->nchar - 2 * qu->trim); 
  else          for (i = 0; i < qu->aln->ntax; i++) non_n[i] = quick_count_sequence_non_N (qu->aln->character->string[i] + qu->trim, qu->aln->nchar - 2 * qu->trim); 

  e = new_empfreq_sort_increasing (non_n, qu->aln->ntax, 2); // 2 => int (other options are char, size_t)
  for (i = 0; i < qu->aln->ntax; i++) non_n[i] = e->i[i].idx; // now non-n will have order index (in original char_vector) of best, second-best, etc.
  char_vector_reorder_strings_from_external_order (qu->aln->character, non_n);
  char_vector_reorder_strings_from_external_order (qu->aln->taxlabel, non_n);
  if (non_n) free (non_n);
  return;
}

void
exclude_redundant_query_sequences (query_t qu, int keep_more_resolved)
{ // this function is independent from reorder_query_structure
  if (!qu->consensus) biomcmc_error ("I can only exclude sequences after indices are created");
  int i, j, *valid, n_valid = 0, dist = 0, red1, red2;

  valid = (int*) biomcmc_malloc (qu->aln->ntax * sizeof (int));
  for (i = 0; i < qu->aln->ntax; i++) valid[i] = 1;

  for (i = 0; i < qu->aln->ntax-1; i++) for (j = i + 1; j < qu->aln->ntax; j++) if (valid[i] && valid[j]) {
    if (qu->acgt)            quick_pairwise_score_acgt (qu->aln->character->string[i], qu->aln->character->string[j], qu->n_idx, 1, &dist, qu->idx);
    else quick_pairwise_score_truncated_idx_indelcheck (qu->aln->character->string[i], qu->aln->character->string[j], qu->n_idx, 1, &dist, qu->idx);
    if (!dist) {  // one of i or j is redundant, and by order, seq[j] has more ACGT than seq[i] (or the same, if they are not equivalent in fact
      if (qu->acgt) { // first go over polymorphic sites (smaller set than non-polymorphic)
        red1 = left_is_resolved_right_acgt (qu->aln->character->string[i], qu->aln->character->string[j], qu->n_idx, qu->idx);
        if (red1 > 1) continue; // i and j are _not_ equivalent (one is not a more resolved version of the other, since both have unique SNPs)
        red2 = left_is_resolved_right_acgt (qu->aln->character->string[i], qu->aln->character->string[j], qu->n_idx_c, qu->idx_c);
      }
      else {
        red1 = left_is_resolved_right (qu->aln->character->string[i], qu->aln->character->string[j], qu->n_idx, qu->idx);
        if (red1 > 1) continue; // i and j are _not_ equivalent (one is not a more resolved version of the other, since both have unique SNPs)
        red2 = left_is_resolved_right (qu->aln->character->string[i], qu->aln->character->string[j], qu->n_idx_c, qu->idx_c);
      }
      if (red2 > 1) continue; // i and j are _not_ equivalent (one is not a more resolved version of the other, since both have unique SNPs)
      if ((!red1) && (!red2)) valid[j] = 0; // 1 and 2 are identical, even at the indels and Ns
      red1 += red2; // now it can be {-2,-1} if i more resolved, {1,2} if j more resolved and zero if not equivalent (red1 and red2 have opposite signs)
      if (!red1) continue; // consensus indices and polymorphic indices have complementary SNP info
      if (keep_more_resolved) {
        if(red1 > 0) valid[i] = 0; // red1 > 0 means that seq[i] is less resolved
        else         valid[j] = 0; // red1 < 0 means that seq[i] is more resolved
      }
      else { 
        if(red1 > 0) valid[j] = 0; // red1 > 0 means that seq[i] is less resolved
        else         valid[i] = 0; // red1 < 0 means that seq[i] is more resolved
      }
    }
  }

  for (i = 0; i < qu->aln->ntax; i++) if (valid[i]) valid[n_valid++] = i; // overwrites, since n_valid is slower than i 
  char_vector_reduce_to_valid_strings (qu->aln->character, valid, n_valid);
  char_vector_reduce_to_valid_strings (qu->aln->taxlabel, valid, n_valid);
  qu->aln->ntax = n_valid;
  if (valid) free (valid);
  return;
}

