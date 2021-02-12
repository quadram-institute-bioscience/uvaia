/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */

#include "utils.h" 

KSEQ_INIT(gzFile, gzread)

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_int  *nbest;
  struct arg_int  *nmax;
  struct arg_int  *trim;
  struct arg_dbl  *ambig;
  struct arg_file *ref;
  struct arg_file *out;
  struct arg_int  *threads;
  struct arg_file *fasta;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .nbest   = arg_int0("n","nbest", NULL, "number of best reference sequences per query to show (default=8)"),
    .nmax    = arg_int0("m","nmax", NULL, "max number of best reference sequences when several optimal (default=2 x nbest)"),
    .trim    = arg_int0(NULL,"trim", NULL, "number of sites to trim from both ends (default=0, suggested for sarscov2=230) -- MAY CONTAIN BUGS"),
    .ambig   = arg_dbl0("a","ambiguity", NULL, "maximum allowed ambiguity for sequence to be excluded (default=0.5)"),
    .ref     = arg_file1("r","reference", "[ref.fa(.gz)]", "*aligned* reference sequences"),
    .out     = arg_file0("o","output", "[chosen_refs.fa.gz]", "GZIPPED output reference sequences (default is to not save sequences)"),
    .threads = arg_int0("t","nthreads",NULL, "suggested number of threads (default is to let system decide; I may not honour your suggestion btw)"),
    .fasta   = arg_filen(NULL, NULL, "[query.fa(.gz)]", 1, 1, "*aligned* sequences to search for neighbour references"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.nbest, params.nmax, params.trim, params.ambig, params.ref, params.out, params.threads, params.fasta, params.end};
  params.argtable = argtable; 
  params.nbest->ival[0] = 8; // default values before parsing
  params.nmax->ival[0] = 0;
  params.trim->ival[0] = 0;
  params.ambig->dval[0] = 0.5;
  /* actual parsing: */
  if (arg_nullcheck(params.argtable)) biomcmc_error ("Problem allocating memory for the argtable (command line arguments) structure");
  arg_parse (argc, argv, params.argtable); // returns >0 if errors were found, but this info also on params.end->count
  print_usage (params, argv[0]);
  return params;
}

// trim for sarscov2: 265 29675

void
del_arg_parameters (arg_parameters params)
{
  if (params.help)  free (params.help);
  if (params.version) free (params.version);
  if (params.nbest) free (params.nbest);
  if (params.nmax) free (params.nmax);
  if (params.trim) free (params.trim);
  if (params.ambig) free (params.ambig);
  if (params.ref)   free (params.ref);
  if (params.fasta) free (params.fasta);
  if (params.out)   free (params.out);
  if (params.threads) free (params.threads);
  if (params.end)   free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }
  if (!params.end->count && (!params.help->count)) return;

  if (params.end->count && (!params.help->count)) {  // params.end holds error messages
    biomcmc_fprintf_colour (stdout, 0,1, "Error when reading arguments from command line:\n", NULL);
    arg_print_errors(stdout, params.end, basename(progname));
  }

  printf ("%s \n", PACKAGE_STRING);
  printf ("Search query sequences against reference ones to describe closest ones\n");
  printf ("This program relies on aligned sequences (queries and references must have same size)\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("Assumes relatively close sequences (e.g. in the SARS-CoV-2 context)\n");
    printf ("Outputs a table with closest sequences to stdout, please redirect as appropriate\n");
  }
  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, *idx = NULL, n_idx = 0;
  int64_t time0[2]; 
  double t_secs = 0.;
  kseq_t *seq;
  gzFile fp;
  size_t trim;
  char_vector cv_seq, cv_name;

  biomcmc_get_time (time0); 
  arg_parameters params = get_parameters_from_argv (argc, argv);

  if (params.nbest->ival[0] < 1) params.nbest->ival[0] = 1;
  if (params.nmax->ival[0] < params.nbest->ival[0]) params.nmax->ival[0] = 2 * params.nbest->ival[0];
  if (params.ambig->dval[0] < 0.001) params.ambig->dval[0] = 0.001;
  if (params.ambig->dval[0] > 1.)    params.ambig->dval[0] = 1.;

#ifdef _OPENMP
  if (params.threads->count) {
    int max_threads = omp_get_max_threads ();
    if (params.threads->ival[0] < 1) params.threads->ival[0] = 1; 
    if (params.threads->ival[0] > max_threads) params.threads->ival[0] = max_threads; 
    omp_set_num_threads (params.threads->ival[0]);
  }
#endif

  /* 1. read reference sequences into char_vectors (ref genome wll be on seq->seq.s and name on seq->name.s) */
  cv_seq  = new_char_vector (1);
  cv_name = new_char_vector (1);

  fp = gzopen((char*) params.ref->filename[0], "r");
  seq = kseq_init(fp);
  while ((i = kseq_read(seq)) >= 0) {
    if (!cv_seq->next_avail) { // first sequence; used here only to fix trim at given value
      if (params.trim->ival[0] <= 0) trim = 0; // lower bound is zero (default value)
      else trim = (size_t) params.trim->ival[0];
      if (trim > seq->seq.l/3) trim = seq->seq.l/3; // upper bound is 1/3 of genome 
    }
    add_reference_genome_to_char_vectors (seq->name.s, seq->seq.s + trim, seq->seq.l - 2 * trim, cv_seq, cv_name, params.ambig->dval[0]);
  }

  gzclose(fp);
  kseq_destroy(seq);
  fprintf (stderr, "finished reading references in %lf secs\n", biomcmc_update_elapsed_time (time0)); fflush(stderr);

  if (cv_seq->nstrings < 1) biomcmc_error ("No valid reference sequences found. Please check file %s.", params.ref->filename[0]);

  /* 2. read each query sequence and align against reference */
  fp = gzopen((char*) params.fasta->filename[0], "r");
  seq = kseq_init(fp); 

  print_score_header ();
  while ((i = kseq_read(seq)) >= 0) 
    t_secs += query_genome_against_char_vectors (seq->name.s, seq->seq.s + trim, seq->seq.l - 2 * trim, cv_seq, cv_name, params.nbest->ival[0], 
                                                 params.nmax->ival[0], &idx, &n_idx, params.ambig->dval[0]);

  fprintf (stderr, "finished search in %lf secs (%lf secs within loop)\n", biomcmc_update_elapsed_time (time0), t_secs); fflush(stderr);

  if (params.out->count) save_sequences (params.out->filename[0], idx, n_idx, cv_seq, cv_name);
  fprintf (stderr, "File saved in %lf secs\n", biomcmc_update_elapsed_time (time0)); fflush(stderr);

  /* everybody is free to feel good */
  kseq_destroy(seq);
  if (idx)   free (idx);
  del_char_vector (cv_seq);
  del_char_vector (cv_name);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}
