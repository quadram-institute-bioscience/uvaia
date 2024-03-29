/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 */


// This is the LEGACY version of uvaia (i.e. used in first publications)
#include "utils.h" 

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_int  *nbest;
  struct arg_int  *nmax;
  struct arg_int  *trim;
  struct arg_dbl  *ambig_r;
  struct arg_dbl  *ambig_q;
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
    .trim    = arg_int0(NULL,"trim", NULL, "number of sites to trim from both ends (default=0, suggested for sarscov2=230)"),
    .ambig_r = arg_dbl0("A","ref_ambiguity",   NULL, "maximum allowed ambiguity for REFERENCE sequence to be excluded (default=0.5)"),
    .ambig_q = arg_dbl0("a","query_ambiguity", NULL, "maximum allowed ambiguity for QUERY sequence to be excluded (default=0.5)"),
    .ref     = arg_file1("r","reference", "[ref.fa(.gz,.xz)]", "*aligned* reference sequences (raw or compressed with gzip, xz, etc.)"),
    .out     = arg_file0("o","output", "[chosen_refs.fa.xz]", "XZIPPED (LZMA) output reference sequences (default is to not save sequences)"),
    .threads = arg_int0("t","nthreads",NULL, "suggested number of threads (default is to let system decide; I may not honour your suggestion btw)"),
    .fasta   = arg_filen(NULL, NULL, "[query.fa(.gz,.xz)]", 1, 1, "*aligned* sequences to search for neighbour references"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.nbest, params.nmax, params.trim, params.ambig_r, params.ambig_q, params.ref, params.out, params.threads, params.fasta, params.end};
  params.argtable = argtable; 
  params.nbest->ival[0] = 8; // default values before parsing
  params.nmax->ival[0] = 0;
  params.trim->ival[0] = 0;
  params.ambig_r->dval[0] = 0.5;
  params.ambig_q->dval[0] = 0.5;
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
  if (params.ambig_r) free (params.ambig_r);
  if (params.ambig_q) free (params.ambig_q);
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
    printf ("Assumes relatively close sequences (e.g. in the SARS-CoV-2 context). ");
    printf ("Outputs a table with closest sequences to stdout, please redirect as appropriate\n");
    printf ("The ambiguities are the highest proportion of Ns allowed; It assumes that the non-Ns are at least 90%% ACGT.\n");

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
  size_t trim;
  //char_vector cv_seq, cv_name;
  alignment refaln, query;

  biomcmc_get_time (time0); 
  arg_parameters params = get_parameters_from_argv (argc, argv);

  if (params.nbest->ival[0] < 1) params.nbest->ival[0] = 1;
  if (params.nmax->ival[0] < params.nbest->ival[0]) params.nmax->ival[0] = 2 * params.nbest->ival[0];
  if (params.ambig_r->dval[0] < 0.001) params.ambig_r->dval[0] = 0.001;
  if (params.ambig_r->dval[0] > 1.)    params.ambig_r->dval[0] = 1.;
  if (params.ambig_q->dval[0] < 0.001) params.ambig_q->dval[0] = 0.001;
  if (params.ambig_q->dval[0] > 1.)    params.ambig_q->dval[0] = 1.;
  
  fprintf (stderr, "Legacy program: %s (historical value only); package: %s\n", basename(argv[0]), PACKAGE_STRING);

#ifdef _OPENMP
  if (params.threads->count) {
    int max_threads = omp_get_max_threads ();
    if (params.threads->ival[0] < 1) params.threads->ival[0] = 1; 
    if (params.threads->ival[0] > max_threads) params.threads->ival[0] = max_threads; 
    omp_set_num_threads (params.threads->ival[0]);
  }
#endif

  /* 1. read reference sequences into char_vectors  */
  refaln = read_fasta_alignment_from_file ((char*)params.ref->filename[0], 0xf); // 0xf -> neither true or false, but gets only bare alignment info

  if (params.trim->ival[0] <= 0) trim = 0; // lower bound is zero (default value)
  else trim = (size_t) params.trim->ival[0];
  if (trim > refaln->nchar / 2.1) trim = refaln->nchar / 2.1; // if we trim more than 1/2 the genome there's nothing left

  fprintf (stderr, "Finished reading %d reference references in %lf secs; now will exclude low quality sequences\n", refaln->ntax, biomcmc_update_elapsed_time (time0)); fflush(stderr);
  uvaia_keep_only_valid_sequences (refaln, params.ambig_r->dval[0], true); // true -> check if seqs are aligned, exiting o.w.
  fprintf (stderr, "Final reference database composed of  %d valid references, spent %lf secs\n", refaln->ntax, biomcmc_update_elapsed_time (time0)); fflush(stderr);
  if (refaln->ntax < 1) biomcmc_error ("No valid reference sequences found. Please check file %s.", params.ref->filename[0]);

  /* 2. read each query sequence and align against reference */
  query = read_fasta_alignment_from_file ((char*)params.fasta->filename[0], 0xf); // 0xf -> neither true or false, but gets only bare alignment info
  fprintf (stderr, "Finished reading %d query references in %lf secs; now will exclude low quality sequences\n", query->ntax, biomcmc_update_elapsed_time (time0)); fflush(stderr);
  uvaia_keep_only_valid_sequences (query, params.ambig_q->dval[0], false); // false -> do not check if seqs are aligned, will be done below w/ reference
  fprintf (stderr, "Final query database composed of  %d valid references, spent %lf secs\n", query->ntax, biomcmc_update_elapsed_time (time0)); fflush(stderr);

  print_score_header ();
   
  for (i = 0; i < query->ntax; i++) 
    t_secs += query_genome_against_char_vectors (query->taxlabel->string[i], query->character->string[i], query->character->nchars[i], refaln->character, refaln->taxlabel, 
                                                 params.nbest->ival[0], params.nmax->ival[0], &idx, &n_idx, trim);

  fprintf (stderr, "finished search in %lf secs (%lf secs within loop)\n", biomcmc_update_elapsed_time (time0), t_secs); fflush(stderr);

  if (params.out->count) {
    save_sequences (params.out->filename[0], idx, n_idx, refaln->character, refaln->taxlabel);
    fprintf (stderr, "File saved in %lf secs\n", biomcmc_update_elapsed_time (time0)); fflush(stderr);
  }

  /* everybody is free to feel good */
  if (idx)   free (idx);
  del_alignment (refaln);
  del_alignment (query);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}
