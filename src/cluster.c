#include <biomcmc.h>

#include "utils.h" 
#include "fastaseq.h"

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_int  *dist;
  struct arg_file *ref;
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
    .dist    = arg_int0("d","distance", NULL, "seqs with less than this SNP differences will be merged (default=1)"),
    .ref     = arg_filen("r","reference", "<ref.fa(.gz,.xz)>", 0, 1, "reference sequence (medoids are furthest from it)"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa(.gz,.xz)>", 1, 1024, "alignments to merge"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.ref, params.fasta, params.end};
  params.argtable = argtable; 
  params.dist->ival[0] = 1;
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
  if (params.dist)  free (params.dist);
  if (params.ref)   free (params.ref);
  if (params.fasta) free (params.fasta);
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
  printf ("Cluster and dedups alignments\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("One-pass clustering equiv to canopy clustering with single, tight distance.");
  }

  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, j;
  int64_t time0[2]; 
  readfasta_t rfas;

  biomcmc_get_time (time0); 
  arg_parameters params = get_parameters_from_argv (argc, argv);
  
  for (j = 0; j < params.fasta->count; j++) {
    rfas = new_readfasta (params.fasta->filename[j]);

    while ((i = readfasta_next (rfas)) >= 0) { // one query per iteration
      printf ("> %s [%lu]\n%s\n",rfas->name, rfas->seqlength, rfas->seq);
    }

    del_readfasta (rfas);
    fprintf (stderr, "Finished reading file %s in %lf secs\n", params.fasta->filename[j], biomcmc_update_elapsed_time (time0)); fflush(stderr);
  }


  /* everybody is free to feel good */
  del_arg_parameters (params);
  return EXIT_SUCCESS;
}
