/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "fastaseq.h"

fastaseq_t
new_fastaseq (void)
{
  fastaseq_t fs = (fastaseq_t) biomcmc_malloc (sizeof (struct fastaseq_struct));
  fs->n_nn = fs->score = 0;
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
  if (fs->name) free (fs->name);
  if (fs->seq)  free (fs->seq);
  free (fs);
  return;
}

void
update_fasta_seq (fastaseq_t to, char **seq, char **name, size_t nchars, int score)
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
  to->score = score;
  to->nchars = nchars; // we assume checks were done outside this function; we accept new sequence length
}

cluster_t
new_cluster (void)
{
  cluster_t clust = (cluster_t) biomcmc_malloc (sizeof (struct cluster_struct));
  clust->fs = NULL;
  clust->n_fs = 0;
  return clust;
}

void
del_cluster (cluster_t clus)
{
  if (!clust) return;
  int i;
  for (i = clust->n_fs-1; i >= 0; i--) del_fastaseq (clust->fs[i]);
  free (clust);
  return;
}

void
add_seq_to_cluster (cluster_t clust, int idx, char **seq, char **name, size_t nchars, int score)
{
  int i;
  if (idx > clust->n_fs) { // new cluster, no similar sequences in cluster set
    idx = clust->n_fs;
    clust->fs = (fastaseq_t*) biomcmc_realloc ((fastaseq_t*) clust->fs, (++clust->n_fs) * sizeof (fastaseq_t));
    clust->fs[idx] = new_fastaseq ();
    update_fasta_seq (clust->fs[idx], *seq, *name, nchars, score);  // STOPHERE
  }

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
  free (rfas);
}

