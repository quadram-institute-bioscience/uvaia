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
  fs->nn = fs->name = fs->seq = NULL;
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
update_fasta_seq (fastaseq_t to, fastaseq_t from)
{
  // 1. replace sequence info 
  if (to->seq) free (to->seq);
  to->seq = from->seq;
  from->seq = NULL;
  // 2. update list of neighbours
  if (to->name) {
    to->nn = (char**) biomcmc_realloc ((char**) to->nn, (to->n_nn + 1) * sizeof (char*));
    to->nn[to->n_nn++] = to->name;
  }
  // 3. replace name
  to->name = from->name;
  from->name = NULL;
  // 4. update socre and seq length info
  to->score = from->score;
  to->nchars = from->nchars;
  del_fastaseq (from);
}

readfasta_t
new_readfasta (char *seqfilename)
{
  readfasta_t rfas = (readfasta_t) biomcmc_malloc (sizeof (struct readfasta_struct));
  seqfile = biomcmc_open_compress (seqfilename, "r");
  rfas->next_name = NULL;
  rfas->fs = new_fastaseq ();
  rfas->line = NULL, rfas->line_read = NULL;
  rfas->linelength = 0;
}

// STOPHERE
int
readfasta_next (readfasta_t rfas)
{
  char *delim = NULL;
  /* start reading file */
  while (biomcmc_getline_compress (&line_read, &linelength, seqfile) != -1) {
    line = line_read; /* the variable *line_read should point always to the same value (no line++ or alike) */
    if (nonempty_fasta_line (line)) { /* each line can be either the sequence or its name, on a strict order */
      /* sequence description (in FASTA jargon); the sequence name */
      if ((delim = strchr (line, '>'))) char_vector_add_string (taxlabel, ++delim);
      /* the sequence itself, which may span several lines */
      else {
        line = remove_space_from_string (line);
        line = uppercase_string (line);
        char_vector_append_string_big_at_position (character, line, taxlabel->next_avail-1); // counter from taxlabel NOT character
      }
    }
  }
  //fclose (seqfile);
  biomcmc_close_compress (seqfile);
  if (line_read) free (line_read);
  char_vector_finalise_big (character);
  align = new_alignment_from_taxlabel_and_character_vectors (taxlabel, character, seqfilename, compact_patterns);
  return align;
}
