#include <biomcmc.h>

#include "fastaseq.h"
#include "min_heap.h"

typedef struct queue_struct* queue_t;

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_lit  *acgt;
  struct arg_lit  *hires;
  struct arg_int  *nbest;
  struct arg_int  *nmax;
  struct arg_int  *trim;
  struct arg_dbl  *ambig_q;
  struct arg_dbl  *ambig_r;
  struct arg_int  *pool;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_file *out;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

struct queue_struct
{
  int max_incompatible, n_query, n_ref, *res; 
  heap_t *heap;
  size_t trim;
  bool *is_best;
  char **seq, **name; // should have a mindist for query seqs, for each element of queue to avoid race conditions
};

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

char* get_outfile_prefix (const char *prefix, size_t *length);

queue_t new_queue (int n_ref, int n_query, int heap_size, size_t trim);
void del_queue (queue_t cq);
void save_sequence_to_compress_stream (file_compress_t out, char *seq, size_t seq_l, char *name, size_t name_l);
void queue_distance_to_consensus (query_t query, queue_t cq, int ir);
void queue_update_min_heaps (query_t query, int iq, queue_t cq, int ir);
void save_distance_table (queue_t cq, query_t query, char *filename);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .acgt    = arg_litn("x","acgt",0, 1, "considers only ACGT sites (i.e. unambiguous SNP differences) in query sequences (experimental)"),
    .hires   = arg_litn("k","keep_resolved",0, 1, "keep more resolved and exclude redundant query seqs (default is to keep all)"),
    .nbest   = arg_int0("n","nbest", NULL, "number of best reference sequences per query to store (default=100)"),
    .trim    = arg_int0("t","trim", NULL, "number of sites to trim from both ends (default=0, suggested for sarscov2=230)"),
    .ambig_q = arg_dbl0("a","query_ambiguity", NULL, "maximum allowed ambiguity for QUERY sequence to be excluded (default=0.5)"),
    .ambig_r = arg_dbl0("A","ref_ambiguity", NULL, "maximum allowed ambiguity for REFERENCE sequence to be excluded (default=0.5)"),
    .pool    = arg_int0("p","pool", NULL, "Pool size, i.e. how many reference seqs are queued to be processed in parallel (larger than number of threads, defaults to 64 per thread)"),
    .ref     = arg_filen("r","reference", "<ref.fa(.gz,.xz)>", 1, 1024, "aligned reference sequences (can be several files)"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa(.gz,.xz)>", 1, 1, "aligned query sequences"),
    .out     = arg_file0("o","output", "<without suffix>", "prefix of xzipped output alignment and table with nearest neighbour sequences"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.acgt, params.hires, params.nbest, params.trim, 
    params.ambig_r, params.ambig_q, params.pool, params.ref, params.fasta, params.out, params.end};
  params.argtable = argtable; 
  params.nbest->ival[0] = 100;
  params.trim->ival[0] = 0;
  params.ambig_r->dval[0] = 0.5;
  params.ambig_q->dval[0] = 0.5;
#ifdef _OPENMP
  params.pool->ival[0] = 64 * omp_get_max_threads (); // default is to have quite a few
#else
  params.pool->ival[0] = 64;
#endif 

  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  arg_parse (argc, argv, params.argtable); // returns >0 if errors were found, but this info also on params.end->count
  print_usage (params, argv[0]);
  return params;
}

void
del_arg_parameters (arg_parameters params)
{
  if (params.help)  free (params.help);
  if (params.version) free (params.version);
  if (params.acgt)   free (params.acgt);
  if (params.hires)  free (params.hires);
  if (params.nbest)  free (params.nbest);
  if (params.trim)   free (params.trim);
  if (params.ambig_q) free (params.ambig_q);
  if (params.ambig_r) free (params.ambig_r);
  if (params.pool) free (params.pool);
  if (params.ref)   free (params.ref);
  if (params.fasta) free (params.fasta);
  if (params.out)   free (params.out);
  if (params.end)   free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }
  if (!params.end->count && (!params.help->count)) return; // regular run

  if (params.end->count && (!params.help->count)) {  // params.end holds error messages
    biomcmc_fprintf_colour (stdout, 0,1, "Error when reading arguments from command line:\n", NULL);
    arg_print_errors(stdout, params.end, basename(progname));
  }

  printf ("%s \n", PACKAGE_STRING);
  printf ("For every query sequence, finds closest neighbours in reference alignment. \n");
  printf ("Notice that this software is multithreaded (and currently there is no control over number of threads)\n"); 
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("Notice that poorly-resolved query sequences (with excess of Ns or indels) have more close neighbours, since only non-N sites are considered). ");
    printf ("You may consider excluding query sequences by chosing a high ");
    printf ("`--query_ambiguity` value especially if they are within the same cluster as another, more resolved query sequence.\n\n");
    printf ("With `--keep_resolved`, if two sequences are equivalent ---in that one is more resolved than another---, only the higher-resolution version is considered. \n\n ");
    printf ("e.g.\nseq1 AAAGCG-C  : higher resolution\nseq2 AA--CG-C  : lower resolution, but distance=0 to seq1 since no SNPs disagree\nseq3 AAA-CGAC  : also no disagreements (d=0) to seq1\n");
    printf ("seq4 AAA-CGTC  : also distance=0 to seq1\n\n");
    printf ("Notice that both seq3 and seq4 have info not available on seq1 and thus are both kept. However seq2 is equivalent to seq1. \n\n");
  }

  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int j, c, n_invalid = 0, non_n_ref, count = 0, n_clust = 256, n_output = 0, print_interval = 50000;
  bool end_of_file = false;
  int64_t time0[2], time1[2];
  double elapsed = 0.;
  size_t outlength = 0;
  char *outfilename = NULL;
  query_t query;  // all query sequences in memory 
  readfasta_t rfas; // single fasta sequence (for queue of ref seqs)
  queue_t cq;       // queue of ref sequences
  file_compress_t outstream;

  biomcmc_get_time (time0); 
  arg_parameters params = get_parameters_from_argv (argc, argv);
  if (params.ambig_q->dval[0] < 0.001) params.ambig_q->dval[0] = 0.001;
  if (params.ambig_q->dval[0] > 1.)    params.ambig_q->dval[0] = 1.;
  if (params.ambig_r->dval[0] < 0.001) params.ambig_r->dval[0] = 0.001;
  if (params.ambig_r->dval[0] > 1.)    params.ambig_r->dval[0] = 1.;
  if (params.nbest->ival[0] < 1)       params.nbest->ival[0] = 1;

#ifdef _OPENMP
  n_clust = omp_get_max_threads (); // upper bound may be distinct to whatever number user has chosen
#else
  n_clust = 1; // compiled without openMP support (e.g. --disable-openmp)
  biomcmc_warning ("Program compiled without multithread support");
#endif
  if (params.pool->ival[0] >= n_clust) n_clust = params.pool->ival[0]; // n_clust was temporary for max number of threads
  fprintf (stderr, "Creating a queue of %d sequences; for each query, the %d closest sequences will be stored\n", n_clust, params.nbest->ival[0]); 
  fflush(stderr);
  
  if (params.out->count) outfilename = get_outfile_prefix (params.out->filename[0], &outlength); // outfilename already has "aln.xz" suffix
  else                   outfilename = get_outfile_prefix ("nn_uvaia", &outlength);

  /* 1. read query sequences into char_vectors */
  query = new_query_structure_from_fasta ((char*)params.fasta->filename[0], params.trim->ival[0], 1, params.acgt->count); // 1=radius distance (not used here)

  fprintf (stderr, "Finished reading %d query references in %lf secs;\n", query->aln->ntax, biomcmc_update_elapsed_time (time0)); fflush(stderr);
  uvaia_keep_only_valid_sequences (query->aln, params.ambig_q->dval[0], true); // true -> check if seqs are aligned, exiting o.w.
  fprintf (stderr, "Query database composed of %d valid references, after excluding low quality.\n", query->aln->ntax); fflush(stderr);
  if (query->aln->ntax < 1) biomcmc_error ("No valid reference sequences found. Please check file %s.", params.fasta->filename[0]);
  biomcmc_get_time (time1); // starts counting here (since we'll have a total time0 and a within-loop time1 chronometer 

  /* 1.1 create indices of polymorphic and monomorphic sites, skipping indels and Ns */
  create_query_indices (query);
  reorder_query_structure (query); // from lower to higher resolution (seqs with more Ns first)

  fprintf (stderr, "Query sequences have %d segregating and %d non-segregating sites (used in comparisons) ", query->n_idx, query->n_idx_c); 
  if (query->acgt) fprintf (stderr, "by focusing on ACGT differences only. \n"); 
  else fprintf (stderr, "by focusing on text match (excluding only indels and Ns).\n"); 
  fflush(stderr);

  if (params.hires->count) { // else do nothing, don't try to remove more resolved and keep less resolved
    exclude_redundant_query_sequences (query, params.hires->count);
    fprintf (stderr, "Updated query database now composed of %d valid references, after removing redundant (less resolved) sequences.\n", query->aln->ntax);
    create_query_indices (query);
    fprintf (stderr, "Query sequences have %d segregating and %d non-segregating sites (both of which are used in comparisons), after removing redundancy\n", query->n_idx, query->n_idx_c); 
    fflush(stderr);
  }

  /* 1.2 several ref sequences with char pointers etc, unrelated to number of threads (just the batch of seqs read at once) */
  cq = new_queue (n_clust, query->aln->ntax, params.nbest->ival[0], query->trim); // query has corrected trim and dist
  /* 1.3 open outfile, trying in order xz, bz, gz, and raw */
  outstream = biomcmc_open_compress (outfilename, "w"); 

  non_n_ref = (int)(query->aln->nchar * params.ambig_r->dval[0]);

  /* 2. read alignment files (can be several) and fill pool of cluster queues */
  fprintf (stderr, "\nNext step is main comparison, which may take a while\n\n"); 
  count = n_invalid = 0;

  for (j = 0; j < params.ref->count; j++) {
    rfas = new_readfasta (params.ref->filename[j]);
    end_of_file = false;

    while (!end_of_file) {
#pragma omp single
      for (c = 0; c < cq->n_ref; c++) { // reads n_clust sequences
        if (readfasta_next (rfas) >= 0) {
          count++;
          if (quick_count_sequence_non_N (rfas->seq, rfas->seqlength) < non_n_ref) {
            if (rfas->seq) free (rfas->seq);
            if (rfas->name) free (rfas->name);
            rfas->seq = rfas->name = cq->seq[c] = cq->name[c] = NULL;
            c--; n_invalid++; // try again this vector element
            continue; // skip to next (c<n_clust) for iteration
          }
          cq->seq[c] = rfas->seq; rfas->seq = NULL;
          cq->name[c] = rfas->name; rfas->name = NULL;
          if (rfas->seqlength != (size_t) query->aln->nchar) {
            biomcmc_warning ("Reference sequence '%s' has %d sites but query sequences have %d sites\n", cq->name[c], rfas->seqlength, query->aln->nchar);
            del_queue (cq);
            del_query_structure (query);
            del_readfasta (rfas);
            biomcmc_error ("all sequences must be aligned");
          }
        }
        else {
          cq->seq[c] = cq->name[c] = NULL;
          end_of_file = true;
        }
      } // single thread for() loop

      // non-threaded
      for (c = 0; c < cq->n_ref * cq->n_query ; c++) cq->is_best[c] = 0; // reset indicators
      cq->max_incompatible = cq->heap[0]->max_incompatible; // we could use a min_heap here, but seems overkill
      for (c = 1; c < cq->n_query; c++) if (cq->max_incompatible < cq->heap[c]->max_incompatible) cq->max_incompatible = cq->heap[c]->max_incompatible;

#pragma omp parallel for shared(c, cq, query) 
      for (c = 0; c < cq->n_ref; c++) { if (cq->seq[c]) queue_distance_to_consensus (query, cq, c); }

#pragma omp parallel for shared(j, query, cq) private (c)
      for (j = 0; j < cq->n_query; j++) {   
        for (c = 0; c < cq->n_ref; c++) if (cq->seq[c]) queue_update_min_heaps (query, j, cq, c);
      }
// Currently, saving seqs is a seq dump. we could alternatively check each ref to see if it's _current_best_match_ for each query, but we migth miss e.g. 
// consecutive sequences that would all belong to best set. Thus better to save more than necessary (any seq that at
// some point was included in the best set, even if ultimately they look far)
#pragma omp parallel for shared(c, cq) private (j)
      for (c = 0; c < cq->n_ref; c++) if (cq->seq[c]) { 
        for (j = 1; j < cq->n_query; j++) cq->is_best[cq->n_query * c + 0] |=  cq->is_best[cq->n_query * c + j]; // onedim [n_ref][n_query]
      }

#pragma omp single // lzma is multithreaded so we must make sure only one thread here
      for (c = 0; c < cq->n_ref; c++) if (cq->seq[c]) { // last round may have fewer sequences than n_ref
        if (cq->is_best[cq->n_query * c] > 0) { // is new best for some of the query seqs
          n_output++;
          save_sequence_to_compress_stream (outstream, cq->seq[c], (size_t) query->aln->nchar, cq->name[c], strlen (cq->name[c]));
        }
      }
      for (c = 0; c < cq->n_ref; c++) { // not used anymore, must be freed by hand here o.w. we have dangling pointers
        if (cq->seq[c]) free (cq->seq[c]);
        if (cq->name[c]) free (cq->name[c]);
        cq->seq[c] = cq->name[c] = NULL;
      }

      if ((count >= print_interval) && ((count % print_interval) < cq->n_ref)) {
        elapsed = biomcmc_update_elapsed_time (time1); 
        fprintf (stderr, "%d sequences analysed in total, %d saved, %d rejected due to high ambiguity; %.3lf secs passed since\n", count, n_output, n_invalid, elapsed); 
        fflush(stderr);
      }

    } // while not end of file

    del_readfasta (rfas);
    fprintf (stderr, "Finished reading file %s in %.3lf secs;\n", params.ref->filename[j], biomcmc_update_elapsed_time (time0));
    fprintf (stderr, "Total of %d sequences read; %d saved sequences include closest neighbours and intermediate, %d too ambiguous (excluded)\n", count, n_output, n_invalid); 
    fflush (stderr);
  }  // for fasta file
  
  fprintf (stderr, "Saved %d sequences to file %s\n", n_output, outstream->filename); 
  biomcmc_close_compress (outstream);

  strcpy (outfilename + outlength, ".csv.xz");
  save_distance_table (cq, query, outfilename);
  fprintf (stderr, "Saved distance table to file %s\n", outfilename); 

  /* everybody is free to feel good */
  del_query_structure (query);
  del_arg_parameters (params);
  del_queue (cq);
  if (outfilename) free (outfilename);
  return EXIT_SUCCESS;
}

char*
get_outfile_prefix (const char *prefix, size_t *length)
{
  char *filename = NULL;
  *length = strlen (prefix);
  filename = (char*) biomcmc_malloc (*length + 8); // suffixes added later
  strncpy (filename, prefix, *length);
  filename[*length] = '\0';
  strcpy (filename + *length, ".aln.xz");  
  return filename;
}
  
queue_t
new_queue (int n_ref, int n_query, int heap_size, size_t trim)
{
  queue_t cq = (queue_t) biomcmc_malloc (sizeof (struct queue_struct));
  int c;
  cq->n_ref = n_ref;
  cq->n_query = n_query;
  cq->trim = trim;
  cq->max_incompatible = 0xffffff;

  cq->res    = (int*)   biomcmc_malloc (4 * n_ref * sizeof (int)); // each ref seq against query->consensus
  cq->seq    = (char**) biomcmc_malloc (n_ref * sizeof (char*)); // only pointer, actual seq is allocated by readfasta
  cq->name   = (char**) biomcmc_malloc (n_ref * sizeof (char*)); // only pointer '' '' ''
  cq->is_best = (bool*) biomcmc_malloc (n_ref * n_query * sizeof (bool*));
  for (c = 0; c < n_ref; c++) { cq->seq[c] = cq->name[c] = NULL;}
  for (c = 0; c < n_ref * n_query ; c++) cq->is_best[c] = 0;
  for (c = 0; c < n_ref * 4; c++) cq->res[c] = 0;

  cq->heap = (heap_t*) biomcmc_malloc (n_query * sizeof (heap_t));
  for (c = 0; c < n_query; c++) cq->heap[c] = new_heap_t (heap_size);

  return cq;
}

void
del_queue (queue_t cq)
{
  if (!cq) return;
  int i;
  if (cq->seq) {
    for (i = cq->n_ref-1; i >= 0; i--) if (cq->seq[i]) free (cq->seq[i]); 
    free (cq->seq);
  }
  if (cq->name) {
    for (i=cq->n_ref-1; i >= 0; i--) if (cq->name[i]) free (cq->name[i]); 
    free (cq->name);
  }
  if (cq->heap) {
    for (i=cq->n_query-1; i >= 0; i--) del_heap_t (cq->heap[i]); 
    free (cq->heap);
  }
  if (cq->is_best) free (cq->is_best);
  if (cq->res) free (cq->res);
  free (cq);
  return;
}

void          
save_sequence_to_compress_stream (file_compress_t xz, char *seq, size_t seq_l, char *name, size_t name_l)
{
  int errors = 0;
  if (biomcmc_write_compress (xz, ">") != 1) errors++; 
  if (biomcmc_write_compress (xz, name) != (int) name_l) errors++;
  if (biomcmc_write_compress (xz, "\n") != 1) errors++;
  if (biomcmc_write_compress (xz, seq) != (int) seq_l) errors++;
  if (biomcmc_write_compress (xz, "\n") != 1) errors++;
  if (errors) fprintf (stderr,"File %s may not be correctly compressed, %d error%s occurred when saving sequence %s.\n", xz->filename, errors, ((errors > 1)? "s":""), name);
}

void
queue_distance_to_consensus (query_t query, queue_t cq, int ir)
{
  biomcmc_pairwise_score_matches_truncated_idx (cq->seq[ir], query->consensus, query->n_idx_c, cq->max_incompatible, cq->res + 4 * ir, (int*) query->idx_c);
}

void
queue_update_min_heaps (query_t query, int iq, queue_t cq, int ir)
{ // compare query iq with reference ir
  q_item this;
  int current_incompatible = cq->res[4*ir + 3] - cq->res[4*ir + 2]; // number of mismatches between this ref and query->consensus;
  this.name = NULL;

  if (current_incompatible > cq->heap[iq]->max_incompatible) return; // this reference is already too far from the query
  current_incompatible = cq->heap[iq]->max_incompatible - current_incompatible; // decrease tolerance, since ref is already "current_incompatible" SNPs away from consensus (and thus, query)

  biomcmc_pairwise_score_matches_truncated_idx (cq->seq[ir], query->aln->character->string[iq], query->n_idx, current_incompatible, this.score, (int*) query->idx);
  if ((this.score[3] - this.score[2]) >= current_incompatible) return; // too far from query

  for (int i = 0; i < 4; i++) this.score[i] += cq->res[4*ir + i]; // score (ref) = score(ref,consensus) + score(ref, this query)
  this.name = cq->name[ir]; // just a pointer (heap_insert copies if necessary)

  if (heap_insert (cq->heap[iq], this)) { // ref ir is within min_heap for query iq
    cq->is_best[cq->n_query * ir + iq] = 1; // onedim [n_ref][n_query]
    if (cq->heap[iq]->n == cq->heap[iq]->heap_size) // heap is full (o.w. the first element might already be best and we'd never fill it)
      cq->heap[iq]->max_incompatible = cq->heap[iq]->seq[1].score[3] - cq->heap[iq]->seq[1].score[0] + 1; // seq[1] contains the worse distance (amongst the best)
  }
} 

void
save_distance_table (queue_t cq, query_t query, char *filename)
{
  file_compress_t xz = biomcmc_open_compress (filename, "w");
  char *buffer = NULL;
  int i, j, k, errors=0;
  size_t buf_len = 128;

  buffer = (char*) biomcmc_malloc (buf_len * sizeof (char));
  sprintf (buffer, "query,reference,rank,ACGT_matches,text_matches,partial_matches,valid_sites\n");
  if (biomcmc_write_compress (xz, buffer) != (int) strlen(buffer)) biomcmc_warning ("problem saving header of compressed file %s;", xz->filename);

  for (i = 0; i < cq->n_query; i++) {
    heap_finalise_heap_qsort (cq->heap[i]); // sort 
    for (j = 0; j < cq->heap[i]->heap_size; j++) {
      buf_len = strlen (cq->heap[i]->seq[j].name) + query->aln->taxlabel->nchars[i] + 64;
      buffer = (char*) biomcmc_realloc ((char*) buffer, buf_len * sizeof (char));
      buffer[0] = '\0';
      sprintf (buffer, "%s,%s,%d", query->aln->taxlabel->string[i], cq->heap[i]->seq[j].name,j+1);
      for (k=0; k < 4; k++) sprintf (buffer + strlen(buffer), ",%d", cq->heap[i]->seq[j].score[k]); // +strlen() to go to last position, o.w. overwrites
      sprintf (buffer + strlen(buffer), "\n");
      if (biomcmc_write_compress (xz, buffer) != (int) strlen(buffer)) errors++; 
      if (buffer) free (buffer);
      buffer = NULL;
    }
  }
  if (errors) fprintf (stderr,"File %s may not have been correctly compressed, %d error%s occurred.\n", xz->filename, errors, ((errors > 1)? "s":""));

  biomcmc_close_compress (xz);
  if (buffer) free (buffer);
}
