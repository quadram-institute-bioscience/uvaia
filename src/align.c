#include "kseq.h"
#include <biomcmc.h>
#include <gap_affine/affine_wavefront_align.h>

KSEQ_INIT(gzFile, gzread)

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);
char * return_query_aligned (char* pattern, int pattern_length, char* text, int text_length, edit_cigar_t* edit_cigar, mm_allocator_t* mm_allocator);
bool sequence_n_below_threshold (char *seq, int seq_length, double threshold);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .ref     = arg_file1("r","reference", "<ref.fa|ref.fa.gz>", "reference sequence"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa|seqs.fa.gz>", 1, 1, "sequences to align"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.ref, params.fasta, params.end};
  params.argtable = argtable; 
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
  if (params.ref)   free (params.ref);
  if (params.fasta) free (params.fasta);
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
  printf ("Align query sequences against a reference\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("Based on the WFA implementation\n");
  }
  del_arg_parameters (params); exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i;
  clock_t time0, time1;
  kseq_t *seq, *ref;
  gzFile fp;
  char *aln_sequence = NULL;

  mm_allocator_t* const mm_allocator = mm_allocator_new(BUFFER_SIZE_8M);
  affine_penalties_t affine_penalties = {.match = 0, .mismatch = 4, .gap_opening = 6, .gap_extension = 2}; // bwa-mem values
  affine_wavefronts_t *affine_wavefronts;

  time0 = clock ();
  arg_parameters params = get_parameters_from_argv (argc, argv);
  
  /* 1. read reference sequence (ref genome wll be on ref->seq.s and name on ref->name.s) */
  fp = gzopen((char*) params.ref->filename[0], "r");
  ref = kseq_init(fp);
  i = kseq_read(ref);
  gzclose(fp);

  //affine_wavefronts = affine_wavefronts_new_complete(ref->seq.l, 2*ref->seq.l, &affine_penalties, NULL, mm_allocator);
  affine_wavefronts = affine_wavefronts_new_reduced (ref->seq.l, 2*ref->seq.l, &affine_penalties, 128, 512, NULL, mm_allocator);

  /* 2. read each query sequence and align against reference */
  fp = gzopen((char*) params.fasta->filename[0], "r");
  seq = kseq_init(fp); 

  while ((i = kseq_read(seq)) >= 0) { // one query per iteration
    if (sequence_n_below_threshold (seq->seq.s, seq->seq.l, 0.5)) { 
      /* 2.1 reset wfa struct and align each seq->seq against ref->seq */
      affine_wavefronts_clear (affine_wavefronts);
      affine_wavefronts_align (affine_wavefronts, ref->seq.s, ref->seq.l, seq->seq.s, seq->seq.l);
      /* 2.2 get alignment itself (wfa gives only cigar) , excluding insertions relative to reference */
      aln_sequence = return_query_aligned (ref->seq.s, ref->seq.l, seq->seq.s, seq->seq.l, &affine_wavefronts->edit_cigar, mm_allocator);
      printf (">%s\n%s\n", seq->name.s, aln_sequence);
      mm_allocator_free (mm_allocator, aln_sequence); // delete aln_sequence memory, that was allocated within allocator
    }
    else fprintf (stderr, "Sequence %s is too ambiguous\n", seq->name.s);
  }

  time1 = clock (); fprintf (stderr, "finished in  %lf secs\n",  (double)(time1-time0)/(double)(CLOCKS_PER_SEC)); fflush(stderr); 

  /* everybody is free to feel good */
  affine_wavefronts_delete (affine_wavefronts);
  mm_allocator_delete (mm_allocator); // dlete allocator itself 
  kseq_destroy(seq);
  kseq_destroy(ref);
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}

char *
return_query_aligned (char* pattern, int pattern_length, char* text, int text_length, edit_cigar_t* edit_cigar, mm_allocator_t* mm_allocator) 
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

