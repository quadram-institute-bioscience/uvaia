#include <biomcmc.h>

#include "fastaseq.h"

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_int  *dist;
  struct arg_int  *trim;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_file *out;
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
    .dist    = arg_int0("d","distance", NULL, "seqs with this SNP differences or less will be merged (default=1)"),
    .trim    = arg_int0("t","trim", NULL, "number of sites to trim from both ends (default=0, suggested for sarscov2=230)"),
    .ref     = arg_filen("r","reference", "<ref.fa(.gz,.xz)>", 0, 1, "reference sequence (medoids are furthest from it)"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa(.gz,.xz)>", 1, 1024, "alignments to merge"),
    .out     = arg_file0("o","output", "<without suffix>", "prefix of xzipped output alignment and cluster table files"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.dist, params.trim, params.ref, params.fasta, params.out, params.end};
  params.argtable = argtable; 
  params.dist->ival[0] = 1;
  params.trim->ival[0] = 0;
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
  if (params.trim) free (params.trim);
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
  printf ("Cluster and dedups alignments\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("One-pass clustering similar to canopy clustering with single, tight distance.\n");
    printf ("If a reference file is provided, its first sequence is used s.t. medoids are updated when resolve a previous one and are farther from it.\n");
    printf ("Otherwise a medoid is the most resolved (less Ns, with no SNP contradicting previous medoid) in cluster\n\n.");
  }

  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, j, count = 0;
  int64_t time0[2], time1[2];
  size_t trim = 0, outlength = 0;
  char *outfilename;
  readfasta_t rfas;
  cluster_t clust;

  biomcmc_get_time (time0); 
  biomcmc_get_time (time1); 
  arg_parameters params = get_parameters_from_argv (argc, argv);
  if (params.dist->ival[0] < 0) params.dist->ival[0] = 0; 

  /* 1. read reference sequence (if missing, use database just to get length) */
  if (params.ref->count) rfas = new_readfasta (params.ref->filename[0]);
  else rfas = new_readfasta (params.fasta->filename[0]);
  readfasta_next (rfas);
  
  if (params.trim->ival[0] > 0) trim = (size_t) params.trim->ival[0];
  if (trim > rfas->seqlength / 2.1) trim = rfas->seqlength / 2.1; // if we trim more than 1/2 the genome there's nothing left
  if (params.dist->ival[0] > rfas->seqlength / 10) params.dist->ival[0] = rfas->seqlength / 10;
  clust = new_cluster (&(rfas->seq), rfas->seqlength, params.dist->ival[0], trim);
  
  del_readfasta (rfas);

  if (!params.ref->count) { // do not use sequence information 
    if (clust->reference) free (clust->reference);
    clust->reference = NULL;
  }
  
  /* 2. read alignment files (can be several) */
  for (j = 0; j < params.fasta->count; j++) {
    rfas = new_readfasta (params.fasta->filename[j]);

    while ((i = readfasta_next (rfas)) >= 0) { // one query per iteration
      check_seq_against_cluster (clust, &(rfas->seq), &(rfas->name), rfas->seqlength);
      if (!(++count % 1000)) 
        fprintf (stderr, "%9d sequences into %9d clusters, spent %lf secs\n", count, clust->n_fs, biomcmc_update_elapsed_time (time1)); fflush(stderr);
    }

    del_readfasta (rfas);
    fprintf (stderr, "Finished reading file %s in %lf secs\n", params.fasta->filename[j], biomcmc_update_elapsed_time (time0)); fflush(stderr);
  }

  if (params.out->count) {
    outlength = strlen (params.out->filename[0]);
    outfilename = (char*) biomcmc_malloc (outlength + 8); // suffixes added later
    strncpy (outfilename, params.out->filename[0], outlength);
    outfilename[outlength] = '\0';
  }
  else {
    outlength = strlen ("cluster_uvaia");
    outfilename = (char*) biomcmc_malloc (outlength + 8); // suffixes added later
    strncpy (outfilename, "cluster_uvaia", outlength);
    outfilename[outlength] = '\0';
  }
 
  qsort (clust->fs, clust->n_fs, sizeof (fastaseq_t), compare_fastaseq);

  strcpy (outfilename + outlength, ".aln.xz");
  save_cluster_to_xz_file (clust, outfilename);

  strcpy (outfilename + outlength, ".csv.xz");
  save_neighbours_to_xz_file (clust, outfilename);
    
  fprintf (stderr, "Finished sorting clusters and saving files in %lf secs\n", biomcmc_update_elapsed_time (time0)); fflush(stderr);

  /* everybody is free to feel good */
  del_arg_parameters (params);
  del_cluster (clust);
  if (outfilename) free (outfilename);
  return EXIT_SUCCESS;
}
