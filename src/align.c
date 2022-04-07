#include <biomcmc.h>
#include <gap_affine/affine_wavefront_align.h>


typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_dbl  *ambig;
  struct arg_int  *pool;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_int  *threads;
  struct arg_file *out;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);
char * return_query_aligned (int pattern_length, char* text, int text_length, edit_cigar_t* edit_cigar, mm_allocator_t* mm_allocator);
bool sequence_n_below_threshold (char *seq, int seq_length, double threshold);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .ambig   = arg_dbl0("a","ambiguity", NULL, "maximum allowed ambiguity for sequence to be excluded (default=0.5)"),
    .pool    = arg_int0("p","pool", NULL, "Pool size, i.e. how many query sequences are queued to be processed in parallel (larger than number of threads, defaults to 64 per thread)"),
    .ref     = arg_file1("r","reference", "<ref.fa|ref.fa.xz>", "reference sequence in fasta format, possibly compressed with gz, xz, bz2"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa|seqs.fa.gz>", 1, 1, "sequences to align"),
    .threads = arg_int0("t","nthreads",NULL, "suggested number of threads (default is to let system decide; I may not honour your suggestion btw)"),
    .out     = arg_file0("o","output", "<without suffix>", "prefix of xzipped output alignment and table with nearest neighbour sequences"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.pool, params.threads, params.out, params.ambig, params.ref, params.fasta, params.end};
  params.argtable = argtable; 
  params.ambig->dval[0] = 0.5;
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
  if (params.threads) free (params.threads);
  if (params.ambig) free (params.ambig);
  if (params.ref)   free (params.ref);
  if (params.fasta) free (params.fasta);
  if (params.pool)  free (params.pool);
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
  printf ("Align query sequences against a reference\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("Based on the WFA implementation https://github.com/smarco/WFA\nOutput is printed to stdout (you should redirect to a file if needed)\n");
    printf ("\n");
  }

  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i;
  int64_t time0[2], time1[2];
  double result[3];
  char *aln_sequence = NULL;
  size_t outlength = 0;
  char *outfilename = NULL;
  readfasta_t rfas, ref;

  mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  affine_penalties_t affine_penalties = {.match = 0, .mismatch = 4, .gap_opening = 6, .gap_extension = 2}; // bwa-mem values
  affine_wavefronts_t *affine_wavefronts;

  arg_parameters params = get_parameters_from_argv (argc, argv);
  biomcmc_get_time (time0); 

  if (params.ambig->dval[0] < 0.001) params.ambig->dval[0] = 0.001;
  if (params.ambig->dval[0] > 1.)    params.ambig->dval[0] = 1.;

#ifdef _OPENMP
  if (params.threads->count) {
    int max_threads = omp_get_max_threads ();
    if (params.threads->ival[0] < 1) params.threads->ival[0] = max_threads; 
    if (params.threads->ival[0] > max_threads) params.threads->ival[0] = max_threads; 
    omp_set_num_threads (params.threads->ival[0]);
  }
#else
  biomcmc_warning ("Program compiled without multithread support");
#endif
  
  if (params.out->count) 
    outfilename = get_outfile_prefix (params.out->filename[0], &outlength); // outfilename already has "aln.xz" suffix
  else {
    char randname[32];
    randname = sprintf ("uvaia.%" PRIu64,time0[1] & 0xffffff);
    outfilename = get_outfile_prefix (randname, &outlength);
  }
  
  /* 1. read reference sequence (ref->seq, ref->name, and ref->seqlength) */
  ref = new_readfasta (params.ref->filename[0]); 
  if (readfasta_next (ref) < 1) biomcmc_error ("Error reading reference sequence %s", params.ref->filename[0]);
  // STOPHERE 

  //affine_wavefronts = affine_wavefronts_new_complete(ref->seq.l, 2*ref->seq.l, &affine_penalties, NULL, mm_allocator);
  affine_wavefronts = affine_wavefronts_new_reduced (ref->seq.l, 2*ref->seq.l, &affine_penalties, 128, 512, NULL, mm_allocator);

  /* 2. read each query sequence and align against reference */
  fp = gzopen((char*) params.fasta->filename[0], "r");
  seq = kseq_init(fp); 

  while ((i = kseq_read(seq)) >= 0) { // one query per iteration
    biomcmc_count_sequence_acgt (seq->seq.s, seq->seq.l, result); 
    if (result[2] > params.ambig->dval[0]) 
      fprintf (stderr, "Sequence %s has proportion of N etc. (=%lf) above threshold of %lf\n", seq->name.s, result[2], params.ambig->dval[0]);
    else if (result[0] < 1. - 1.1 * params.ambig->dval[0]) 
      fprintf (stderr, "Sequence %s has proportion of ACGT (=%lf) below threshold of %lf\n", seq->name.s, result[0], 1. - 1.1 * params.ambig->dval[0]);
    else {    //if (sequence_n_below_threshold (seq->seq.s, seq->seq.l, params.ambig->dval[0])) { 
      /* 2.1 reset wfa struct and align each seq->seq against ref->seq */
      affine_wavefronts_clear (affine_wavefronts);
      affine_wavefronts_align (affine_wavefronts, ref->seq.s, ref->seq.l, seq->seq.s, seq->seq.l);
      /* 2.2 get alignment itself (wfa gives only cigar) , excluding insertions relative to reference */
      aln_sequence = return_query_aligned (ref->seq.l, seq->seq.s, seq->seq.l, &affine_wavefronts->edit_cigar, mm_allocator);
      printf (">%s\n%s\n", seq->name.s, aln_sequence);
      mm_allocator_free (mm_allocator, aln_sequence); // delete aln_sequence memory, that was allocated within allocator
    }
  }

  time1 = clock (); fprintf (stderr, "finished in  %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); 

  /* everybody is free to feel good */
  affine_wavefronts_delete (affine_wavefronts);
  mm_allocator_delete (mm_allocator); // dlete allocator itself 
  del_arg_parameters (params);
  if (outfilename) free (outfilename);
  del_readfasta (ref);
  del_readfasta (rfas);
  return EXIT_SUCCESS;
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
  /*i=0;
  while (text_pos < text_length) {
    text_alg[alg_pos+i] = text[text_pos++];
    ops_alg[alg_pos+i] = '?';
    ++i;
  } */
  text_alg[alg_pos] = '\0';
  return text_alg;
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

// FIXME: below is copy-paste from nearest.c
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
new_queue (int n_ref, int n_query, int heap_size, size_t trim, int n_sites)
{
  queue_t cq = (queue_t) biomcmc_malloc (sizeof (struct queue_struct));
  int c;
  cq->n_ref = n_ref;
  cq->n_query = n_query;
  cq->trim = trim;
  cq->max_incompatible = n_sites; // can be a lage number, this is more interpretable as output to user

  cq->res    = (int*)   biomcmc_malloc (4 * n_ref * sizeof (int)); // each ref seq against query->consensus; OBS: acgt mode could use only 2 
  cq->non_n  = (int*)   biomcmc_malloc (n_ref * sizeof (int)); // non-n (i.e. ACGT etc.) sites 
  cq->seq    = (char**) biomcmc_malloc (n_ref * sizeof (char*)); // only pointer, actual seq is allocated by readfasta
  cq->name   = (char**) biomcmc_malloc (n_ref * sizeof (char*)); // only pointer '' '' ''
  cq->is_best = (bool*) biomcmc_malloc (n_ref * n_query * sizeof (bool*));
  for (c = 0; c < n_ref; c++) { cq->seq[c] = cq->name[c] = NULL; cq->non_n[c] = 0;}
  for (c = 0; c < n_ref * n_query ; c++) cq->is_best[c] = 0;
  for (c = 0; c < n_ref * 4; c++) cq->res[c] = 0;

  cq->heap = (heap_t*) biomcmc_malloc (n_query * sizeof (heap_t));
  for (c = 0; c < n_query; c++) { cq->heap[c] = new_heap_t (heap_size); cq->heap[c]->max_incompatible = cq->max_incompatible;}

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
  if (cq->non_n) free (cq->non_n);
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
{  // consensus may have "AAAA" where a query has "----" thus we restrict to idx_c which are identical, non-gap in all queries
  if (query->acgt) quick_pairwise_score_acgt_and_valid (cq->seq[ir], query->consensus, query->n_idx_c, cq->max_incompatible, cq->res + 4 * ir, query->idx_c);
  else    biomcmc_pairwise_score_matches_truncated_idx (cq->seq[ir], query->consensus, query->n_idx_c, cq->max_incompatible, cq->res + 4 * ir, query->idx_c);
}

void
queue_update_min_heaps (query_t query, int iq, queue_t cq, int ir)
{
  if (query->acgt) queue_update_min_heaps_acgt (query, iq, cq, ir);
  else             queue_update_min_heaps_full (query, iq, cq, ir);
}

void
queue_update_min_heaps_acgt (query_t query, int iq, queue_t cq, int ir)
{ // compare query iq with reference ir; result[] will store matches (to reject terrible quality sequences with few mismatches) but quick_pairwise returns mismatches
  q_item this;
  int result[4];
  int current_incompatible = cq->res[4*ir]; // number of mismatches between this reference and query->consensus (acgt result[] needs only two scores) 
  this.name = NULL;

  if (current_incompatible >= cq->heap[iq]->max_incompatible) return; // this reference is already too far from the query
  current_incompatible = cq->heap[iq]->max_incompatible - current_incompatible; // decrease tolerance, since ref is already "current_incompatible" SNPs away from consensus (and thus, query)

  quick_pairwise_score_acgt_and_valid (cq->seq[ir], query->aln->character->string[iq], query->n_idx_m, current_incompatible, result, query->idx_m); // some queries have indels
  if (result[0] >= current_incompatible) return; // this reference is already too far from the query
  current_incompatible = current_incompatible - result[0]; // decrease tolerance, and result[0] has mismatches
  result[0] += cq->res[4*ir];     // accumulate consensus scores (cq->res[] has score common to all queries, result[] has scores which may involve gaps)
  result[1] += cq->res[4*ir + 1];

  quick_pairwise_score_acgt_and_valid (cq->seq[ir], query->aln->character->string[iq], query->n_idx, current_incompatible, result + 2, query->idx);
  if (result[2] >= current_incompatible) return; 

  /* we accumulate first 2 scores between this query and consensus (the 2 scores come from quick_pairwise_acgt) */
  for (int i = 0; i < 6; i++) this.score[i] = 0; 
  this.score[0] = result[1] + result[3] - result[0] - result[2]; // result[0] and [2] are mismatches, but we need matches to avoid the low quality trap 
  this.score[1] = result[1] + result[3];  // score (ref) = score(ref,consensus) + score(ref, this query)
  this.score[2] = result[3] - result[2]; // unique matches at this sequence (so refs more distant from consensus query are preferred)
  this.score[3] = cq->non_n[ir]; // tie-breaker is the number of valid sites (i.e. excluding indels and Ns) WARNING: can bias towards labs which impute sites
  this.score[4] = result[0]; // min_heap unlikely to come here, but good to print: number of mismatches to consensus 
  this.score[5] = result[2]; // min_heap unlikely to come here, but good to print: number of mismatches to this query sequence
  this.name = cq->name[ir];  // just a pointer (heap_insert copies if necessary)

  if (heap_insert (cq->heap[iq], this)) { // ref ir is within min_heap for query iq
    cq->is_best[cq->n_query * ir + iq] = 1; // onedim [n_ref][n_query]; this will make this ref seq being saved to output file
    if (cq->heap[iq]->n == cq->heap[iq]->heap_size) // heap is full (o.w. the first element might already be best and we'd never fill it)
      cq->heap[iq]->max_incompatible = cq->heap[iq]->seq[1].score[1] - cq->heap[iq]->seq[1].score[0] + 1; // seq[1] contains the worse distance (amongst the best); we need mismatches again
  }
} 

void
queue_update_min_heaps_full (query_t query, int iq, queue_t cq, int ir)
{ // compare query iq with reference ir
  q_item this;
  int result[8];
  //int current_incompatible = query->n_idx + query->n_idx_c - cq->res[4*ir + 0]; // number of mismatches between this ref and query->consensus;
  int current_incompatible = cq->res[4*ir + 3] - cq->res[4*ir + 0]; // number of mismatches between this ref and query->consensus;
  this.name = NULL;

  if (current_incompatible >= cq->heap[iq]->max_incompatible) return; // this reference is already too far from the query
  current_incompatible = cq->heap[iq]->max_incompatible - current_incompatible; // decrease tolerance, since ref is already "current_incompatible" SNPs away from consensus (and thus, query)

  biomcmc_pairwise_score_matches_truncated_idx (cq->seq[ir], query->aln->character->string[iq], query->n_idx_m, current_incompatible, result, query->idx_m);
  if ((result[3] - result[0]) >= current_incompatible) return; // too far from query (BTW r[3]-r[0] is same equation used by biomcmc_pairwise())
  current_incompatible = current_incompatible - result[3] + result[0]; // decrease tolerance

  biomcmc_pairwise_score_matches_truncated_idx (cq->seq[ir], query->aln->character->string[iq], query->n_idx, current_incompatible, result + 4, query->idx);
  if ((result[7] - result[4]) >= current_incompatible) return; // too far from query 

  /* we accumulate first 4 scores between this query and two sets from consensus (each 4 scores come from biomcmc_pairwise) */
  for (int i = 0; i < 4; i++) this.score[i] = result[i] + result[i+4] + cq->res[4*ir + i]; // score (ref) = score(ref,consensus) + score(ref, this query)
  this.score[4] = result[0]; // 5th score is unique matches at this sequence (so refs more distant from consensus query are preferred)  
  this.score[5] = cq->non_n[ir]; // last (6th) score is the number of valid sites (i.e. excluding indels and Ns) WARNING: can bias towards labs which impute sites
  this.name = cq->name[ir]; // just a pointer (heap_insert copies if necessary)

  if (heap_insert (cq->heap[iq], this)) { // ref ir is within min_heap for query iq
    cq->is_best[cq->n_query * ir + iq] = 1; // onedim [n_ref][n_query]
    if (cq->heap[iq]->n == cq->heap[iq]->heap_size) // heap is full (o.w. the first element might already be best and we'd never fill it)
      //cq->heap[iq]->max_incompatible = query->n_idx + query->n_idx_c - cq->heap[iq]->seq[1].score[0] + 1; // seq[1] contains the worse distance (amongst the best)
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
  if (query->acgt)
    sprintf (buffer, "query,reference,rank,ACGT_matches,valid_ACGT_comparisons,ACGT_matches_unique,valid_ref_sites,dist_consensus,dist_unique\n");
  else 
    sprintf (buffer, "query,reference,rank,ACGT_matches,text_matches,partial_matches,valid_pair_comparisons,ACGT_matches_unique,valid_ref_sites\n");

  if (biomcmc_write_compress (xz, buffer) != (int) strlen(buffer)) biomcmc_warning ("problem saving header of compressed file %s;", xz->filename);

  for (i = 0; i < cq->n_query; i++) {
    heap_finalise_heap_qsort (cq->heap[i]); // sort 
    for (j = 0; j < cq->heap[i]->heap_size; j++) {
      buf_len = strlen (cq->heap[i]->seq[j].name) + query->aln->taxlabel->nchars[i] + 80;
      buffer = (char*) biomcmc_realloc ((char*) buffer, buf_len * sizeof (char));
      buffer[0] = '\0';
      sprintf (buffer, "%s,%s,%d", query->aln->taxlabel->string[i], cq->heap[i]->seq[j].name,j+1);
      for (k=0; k < 6; k++) sprintf (buffer + strlen(buffer), ",%d", cq->heap[i]->seq[j].score[k]); // +strlen() to go to last position, o.w. overwrites
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
