#include "kseq.h"
#include <biomcmc.h>
#include <gap_affine/affine_wavefront_align.h>

KSEQ_INIT(gzFile, gzread)

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_int  *nbest;
  struct arg_int  *nmax;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_file *out;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);
void upper_kseq (kseq_t *seq);

void describe_scores (char *query_name, double *score, char_vector refnames, int nbest, int nmax, int **idx, int *n_idx);
void save_sequences (const char *filename, int *idx, int n_idx, char_vector seq, char_vector name);
bool sequence_n_below_threshold (char *seq, int seq_length, double threshold);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .nbest   = arg_int0("n","nbest", NULL, "number of best reference sequences per query to show (default=8)"),
    .nmax    = arg_int0("m","nmax", NULL, "max number of best reference sequences when several optimal (default=2 x nbest)"),
    .ref     = arg_file1("r","reference", "<ref.fa|ref.fa.gz>", "_aligned_ reference sequences"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa|seqs.fa.gz>", 1, 1, "_aligned_ sequences to search on references"),
    .out     = arg_file0("o","output", "<gzipped fasta>", "output reference sequences"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.nbest, params.nmax, params.ref, params.fasta, params.out, params.end};
  params.argtable = argtable; 
  params.nbest->ival[0] = 8; // default values before parsing
  params.nmax->ival[0] = 0;  // default values before parsing
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
  if (params.nbest) free (params.nbest);
  if (params.nmax) free (params.nmax);
  if (params.ref)   free (params.ref);
  if (params.fasta) free (params.fasta);
  if (params.out)   free (params.out);
  if (params.end)   free (params.end);
}

void 
print_usage (arg_parameters params, char *progname)
{
  if (params.version->count) { printf ("%s\n", PACKAGE_VERSION); del_arg_parameters (params); exit (EXIT_SUCCESS); }
  if (!params.end->count && (!params.help->count)) return;

  if (params.end->count) {  // params.end holds error messages
    arg_print_errors(stdout, params.end, basename(progname));
    printf ("Error when reading arguments from command line\n\n");
  }

  printf ("%s \n", PACKAGE_STRING);
  printf ("Search query sequences against reference ones to describe closest ones\n");
  printf ("This program relies on aligned sequences (queries and references must have same size)\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("Somehow assumes close sequences (it's being developed in the SARS-CoV-2 context)\n");
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, j, *idx = NULL, n_idx = 0;
  double *score = NULL, k2p[3];
  clock_t time0, time1;
  kseq_t *seq;
  gzFile fp;
  char *aln_sequence = NULL;
  char_vector cv_seq, cv_name;

  time0 = clock ();
  arg_parameters params = get_parameters_from_argv (argc, argv);

  if (params.nbest->ival[0] < 1) params.nbest->ival[0] = 1;
  if (params.nmax->ival[0] < params.nbest->ival[0]) params.nmax->ival[0] = 2 * params.nbest->ival[0];

  /* 1. read reference sequences into char_vectors (ref genome wll be on seq->seq.s and name on seq->name.s) */
  cv_seq  = new_char_vector (1);
  cv_name = new_char_vector (1);
  fp = gzopen((char*) params.ref->filename[0], "r");
  seq = kseq_init(fp);
  while ((i = kseq_read(seq)) >= 0) { // one query per iteration
    if (cv_seq->next_avail && (cv_seq->nchars[cv_seq->next_avail-1] != seq->seq.l))
      biomcmc_error ("This program assumes aligned sequences, and sequence %s has length %u while sequence %s has length %u\n",
                     seq->name.s, seq->seq.l, cv_name->string[cv_seq->next_avail-1], cv_seq->nchars[cv_seq->next_avail-1]);
    if (sequence_n_below_threshold (seq->seq.s, seq->seq.l, 0.5)) { 
      upper_kseq (seq);
      char_vector_add_string (cv_seq,  seq->seq.s);
      char_vector_add_string (cv_name, seq->name.s);
    }
    else fprintf (stderr, "Reference sequence %s is too ambiguous\n", seq->name.s);
  }
  gzclose(fp);
  kseq_destroy(seq);
  time1 = clock (); fprintf (stderr, "finished reading references in %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); time0 = time1; 

  score = (double*) biomcmc_malloc (cv_seq->nstrings * sizeof (double));
  for (j = 0; j < cv_seq->nstrings; j++) score[j] = 0.;

  /* 2. read each query sequence and align against reference */
  fp = gzopen((char*) params.fasta->filename[0], "r");
  seq = kseq_init(fp); 

  while ((i = kseq_read(seq)) >= 0) { // one query per iteration
    if (sequence_n_below_threshold (seq->seq.s, seq->seq.l, 0.5)) { 
      upper_kseq (seq);
#ifdef _OPENMP
#pragma omp parallel for shared(cv_seq,cv_name,seq) schedule(dynamic)
#endif
      for (j = 0; j < cv_seq->nstrings; j++) {
        if (cv_seq->nchars[j] != seq->seq.l)
          biomcmc_error ("This program assumes aligned sequences, and query sequence %s has length %u while reference sequence %s has length %u\n",
                         seq->name.s, seq->seq.l, cv_name->string[j], cv_seq->nchars[j]);
        score[j] = biomcmc_pairwise_score_matches (cv_seq->string[j], seq->seq.s, seq->seq.l, k2p);
      }
      describe_scores (seq->name.s, score, cv_name, params.nbest->ival[0], params.nmax->ival[0], &idx, &n_idx);
    }
    else fprintf (stderr, "Query sequence %s is too ambiguous\n", seq->name.s);
  }
  time1 = clock (); fprintf (stderr, "finished search in %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); time0 = time1;  

  if (params.out->count) save_sequences (params.out->filename[0], idx, n_idx, cv_seq, cv_name);
  time1 = clock (); fprintf (stderr, "finished saving in %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); time0 = time1;  

  /* everybody is free to feel good */
  kseq_destroy(seq);
  if (score) free (score);
  if (idx)   free (idx);
  del_char_vector (cv_seq);
  del_char_vector (cv_name);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}

void
upper_kseq (kseq_t *seq)
{
  for (int i = 0; i < seq->seq.l; i++) seq->seq.s[i] = toupper(seq->seq.s[i]);
}

void 
describe_scores (char *query_name, double *score, char_vector refnames, int nbest, int nmax, int **idx, int *n_idx)
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

