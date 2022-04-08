#include <biomcmc.h>
#include "fastaseq.h"
#include <gap_affine/affine_wavefront_align.h>

typedef struct queue_struct* queue_t;

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_lit  *screen;
  struct arg_dbl  *ambig;
  struct arg_int  *pool;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_int  *threads;
  struct arg_file *out;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

struct queue_struct
{
  int n_seq, n_threads;
  size_t aln_length;
  affine_wavefronts_t **affine_wavefronts; // array of pointers to affine_wavefronts_t
  char *refseq, **seq, **aln, **name; 
};

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

char* get_outfile_prefix (const char *prefix, size_t *length);
queue_t new_queue (char **refseq, size_t n_sites, int n_seq, int n_threads);
void del_queue (queue_t cq);
void save_sequence_to_compress_stream (file_compress_t xz, char *seq, size_t seq_l, char *name, size_t name_l);
void align_query (queue_t cq, int c, int thread);
void update_query_aligned (edit_cigar_t* edit_cigar, queue_t qc, int location);


arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .screen  = arg_litn(NULL,"stdout",0, 1, "print alignment to stdout (to redirect/pipe) instead of compress to file; much faster but may generate a big output"),
    .ambig   = arg_dbl0("a","ambiguity", NULL, "maximum allowed ambiguity for sequence to be excluded (default=0.5)"),
    .pool    = arg_int0("p","pool", NULL, "Pool size, i.e. how many query sequences are queued to be aligned in parallel (larger than number of threads, defaults to 64 per thread)"),
    .ref     = arg_file1("r","reference", "<ref.fa|ref.fa.xz>", "reference sequence in fasta format, possibly compressed with gz, xz, bz2"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa|seqs.fa.xz>", 1, 1024, "sequences to align in fasta format, possibly compressed with gz, xz, bz2 (can be multiple files)"),
    .threads = arg_int0("t","nthreads",NULL, "suggested number of threads (default is to let system decide; I may not honour your suggestion btw)"),
    .out     = arg_file0("o","output", "<without suffix>", "prefix of xzipped output alignment and table with nearest neighbour sequences"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.screen, params.pool, params.threads, params.out, params.ambig, params.ref, params.fasta, params.end};
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
    printf ("Based on the WFA implementation https://github.com/smarco/WFA\nOutput is printed to stdout (you should redirect to a file if needed).\n");
    printf ("Since the sequences are assumed to be similar, sequences too short or too big w.r.t. the reference are rejected.\n\n");
    printf ("The reference sequence and the unaligned fasta files can be compressed with gz, xz, bz2. The alignment output will be compressed with xz, ");
    printf ("unless you miss the library or the software was compiled without support. In this case the next available compression is tried.\n\n");
  }

  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, j, c, count = 0, n_output = 0, print_interval = 5000;
  int64_t time0[2], time1[2];
  double result[3], elapsed;
  bool seq_valid = true, end_of_file = false;
  size_t outlength = 0;
  char *outfilename = NULL;
  readfasta_t rfas;
  file_compress_t outstream;
  queue_t cq = NULL;

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

  if (params.screen->count) 
    fprintf (stderr, "Sequences will be shown uncompressed in screen (to redirect to file or pipe into another software).\n");
  else {
    if (params.out->count) 
      outfilename = get_outfile_prefix (params.out->filename[0], &outlength); // outfilename already has "aln.xz" suffix
    else {
      char randname[32];
      sprintf (randname, "uvaia.%" PRIx64,time0[1] & 0xffffff);
      outfilename = get_outfile_prefix (randname, &outlength);
    }
    fprintf (stderr, "Sequences will be compressed (if possible) and saved into file %s.\n", outfilename);
  }
  
  /* 1. read reference sequence (ref->seq, ref->name, and ref->seqlength) */
  rfas = new_readfasta (params.ref->filename[0]); 
  if (readfasta_next (rfas) < 1) biomcmc_error ("Error reading reference sequence %s", params.ref->filename[0]);

  /* 1.1 create pool of unaligned query sequences which are aligned against the reference in parallel */
#ifdef _OPENMP
  i = omp_get_max_threads ();
#else
  i = 1;
#endif
  cq = new_queue (&(rfas->seq), rfas->seqlength, params.pool->ival[0], i);
  rfas->seq = NULL; // pointer copied above, we don't want to free it here
  del_readfasta (rfas); // will recycle rfas variable to read query unaligned sequences
  /* 1.2 open output directory with aligned sequences */
  if (!params.screen->count) outstream = biomcmc_open_compress (outfilename, "w"); 
  
  biomcmc_get_time (time1); 
  for (j = 0; j < params.fasta->count; j++) {
    fprintf (stderr, "Started  reading file %s\n", params.fasta->filename[j]);
    rfas = new_readfasta (params.fasta->filename[j]);
    end_of_file = false;

    while (!end_of_file) {
#pragma omp single
      for (c = 0; c < cq->n_seq; c++) { 
        if (readfasta_next (rfas) >= 0) {
          count++;
          seq_valid = true;

          if (((3 * rfas->seqlength) < (2 * cq->aln_length)) || ((2 * rfas->seqlength) > (3 * cq->aln_length))) {
            fprintf (stderr, "Sequence %s has size too different from reference (%lu vs %lu)\n", rfas->name, rfas->seqlength, cq->aln_length);
            seq_valid = false;
          }
          if (seq_valid) biomcmc_count_sequence_acgt (rfas->seq, rfas->seqlength, result); 
          if (seq_valid && (result[2] > params.ambig->dval[0])) {
            fprintf (stderr, "Sequence %s has proportion of N etc. (=%lf) above threshold of %lf\n", rfas->name, result[2], params.ambig->dval[0]);
            seq_valid = false;
          }
          if (seq_valid && (result[0] < 1. - 1.1 * params.ambig->dval[0])) { 
            fprintf (stderr, "Sequence %s has proportion of ACGT (=%lf) below threshold of %lf\n", rfas->name, result[0], 1. - 1.1 * params.ambig->dval[0]);
            seq_valid = false;
          }

          if (!seq_valid) {
            if (rfas->seq) free (rfas->seq);
            if (rfas->name) free (rfas->name);
            rfas->seq = rfas->name = cq->seq[c] = cq->name[c] = NULL;
            c--; // try again this vector element
            continue; // skip to next for(c < n_seq) iteration
          }

          cq->seq[c] = rfas->seq; rfas->seq = NULL;
          cq->name[c] = rfas->name; rfas->name = NULL;
        } else { // readfasta < 0 means end of file
          cq->seq[c] = cq->name[c] = NULL;
          end_of_file = true;
        }
      } // single thread for() loop

#pragma omp parallel for shared(c, cq) private(i)
      for (c = 0; c < cq->n_seq; c++) if (cq->seq[c]) {
#ifdef _OPENMP
        i = omp_get_thread_num ();
#else
        i = 0;
#endif
        align_query (cq, c, i);
      }

      if (params.screen->count) {
        for (c = 0; c < cq->n_seq; c++) if (cq->seq[c]) printf (">%s\n%s\n", cq->name[c], cq->aln[c]);
      }
      else {
#pragma omp single // lzma is multithreaded so we must make sure only one thread here
        for (c = 0; c < cq->n_seq; c++) if (cq->seq[c]) { // last round may have fewer sequences than n_ref
          n_output++;
          save_sequence_to_compress_stream (outstream, cq->aln[c], (size_t) cq->aln_length, cq->name[c], strlen (cq->name[c]));
        }
      } // else (save to file)
      for (c = 0; c < cq->n_seq; c++) { // not used anymore, must be freed by hand here o.w. we have dangling pointers
        if (cq->seq[c]) free (cq->seq[c]);
        if (cq->name[c]) free (cq->name[c]);
        cq->seq[c] = cq->name[c] = NULL;
      }

      if ((count >= print_interval) && ((count % print_interval) < cq->n_seq)) {
        elapsed = biomcmc_update_elapsed_time (time1); 
        fprintf (stderr, "%d\t sequences read, %d \t aligned. %.3lf secs elapsed.\n", count, n_output, elapsed);
        fflush(stderr);
      }

    } // while not end of file
    del_readfasta (rfas);
    fprintf (stderr, "Finished reading file %s in %.3lf secs;\n", params.fasta->filename[j], biomcmc_update_elapsed_time (time1));
    fflush (stderr);
  }  // for fasta file
  
  if (params.screen->count) {
    fprintf (stderr, "Output %d aligned sequences. Total elapsed time: %.3lf secs\n", n_output, biomcmc_update_elapsed_time (time0)); 
  }
  else {
    biomcmc_close_compress (outstream);
    fprintf (stderr, "Saved %d sequences to file %s\nTotal elapsed time: %.3lf secs\n", n_output, outfilename,  biomcmc_update_elapsed_time (time0)); 
  }

  /* everybody is free to feel good */
  del_arg_parameters (params);
  if (outfilename) free (outfilename);
  del_queue (cq);
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
new_queue (char **refseq, size_t n_sites, int n_seq, int n_threads)
{
  queue_t cq = (queue_t) biomcmc_malloc (sizeof (struct queue_struct));
  int c;

  cq->n_seq = n_seq;
  cq->n_threads = n_threads;
  cq->aln_length = n_sites;
  cq->refseq = *refseq;

  cq->seq    = (char**) biomcmc_malloc (n_seq * sizeof (char*)); // only pointer, actual seq is allocated by readfasta
  cq->name   = (char**) biomcmc_malloc (n_seq * sizeof (char*)); // only pointer '' '' ''
  cq->aln    = (char**) biomcmc_malloc (n_seq * sizeof (char*)); // actually allocated here since all same length
  for (c = 0; c < n_seq; c++) { 
    cq->seq[c] = cq->name[c] = NULL;
    cq->aln[c] = biomcmc_malloc (n_sites * sizeof (char));
  }
  
  cq->affine_wavefronts = (affine_wavefronts_t**) biomcmc_malloc (n_threads * sizeof (affine_wavefronts_t*));
  
  affine_penalties_t affine_penalties = {.match = 0, .mismatch = 4, .gap_opening = 6, .gap_extension = 2}; // bwa-mem values
  for (c = 0; c < n_threads; c++) {
    mm_allocator_t* const mm_allocator = mm_allocator_new (BUFFER_SIZE_8M); // memory allocator from FWA
    cq->affine_wavefronts[c] = affine_wavefronts_new_reduced (n_sites, 3 * n_sites, &affine_penalties, 128, 512, NULL, mm_allocator);
//    mm_allocator_delete (mm_allocator); // delete allocator itself (cannot be deleted since it's pointed to by affine_wavefronts)
  }

  return cq;
}

void
del_queue (queue_t cq)
{
  if (!cq) return;
  int i;
  if (cq->seq) {
    for (i = cq->n_seq-1; i >= 0; i--) if (cq->seq[i]) free (cq->seq[i]); 
    free (cq->seq);
  }
  if (cq->name) {
    for (i=cq->n_seq-1; i >= 0; i--) if (cq->name[i]) free (cq->name[i]); 
    free (cq->name);
  }
  if (cq->aln) {
    for (i=cq->n_seq-1; i >= 0; i--) if (cq->aln[i]) free (cq->aln[i]); 
    free (cq->aln);
  }
  if (cq->affine_wavefronts) {
    for (i=cq->n_threads-1; i >= 0; i--) if (cq->affine_wavefronts[i]) {
      //mm_allocator_delete (cq->affine_wavefronts[i]->mm_allocator);
      affine_wavefronts_delete (cq->affine_wavefronts[i]);
    }
    free (cq->affine_wavefronts);
  }
  if (cq->refseq) free (cq->refseq);
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
align_query (queue_t cq, int c, int thread)
{
  affine_wavefronts_clear (cq->affine_wavefronts[thread]);
  affine_wavefronts_align (cq->affine_wavefronts[thread], cq->refseq, cq->aln_length, cq->seq[c], strlen (cq->seq[c]));
  update_query_aligned (&(cq->affine_wavefronts[thread]->edit_cigar), cq, c);
  return;
}

void
update_query_aligned (edit_cigar_t* edit_cigar, queue_t qc, int location) 
{
  // Parameters
  char* const operations = edit_cigar->operations;
  // Compute alignment buffers
  int i, alg_pos = 0, pattern_pos = 0, text_pos = 0;
  for (i = edit_cigar->begin_offset; i < edit_cigar->end_offset; ++i) switch (operations[i]) {
    case 'M':
    case 'X':
      qc->aln[location][alg_pos++] = qc->seq[location][text_pos++];
      pattern_pos++;
      break;
    case 'I':
      text_pos++;
      break;
    case 'D':
      qc->aln[location][alg_pos++] = '-';
      pattern_pos++;
      break;
    default:
      break;
  }
  qc->aln[location][alg_pos++] = '\0';
}

