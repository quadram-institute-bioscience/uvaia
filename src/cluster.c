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
  int j, c, count = 0, n_threads = 1, print_interval = 500;
  bool end_of_file = false;
  int64_t time0[2], time1[2];
  size_t trim = 0, outlength = 0, *nchars_vec;
  char *outfilename = NULL, *refseq = NULL, **seq_vec, **name_vec;
  readfasta_t rfas;
  cluster_t *clust;

  #ifdef _OPENMP
  n_threads = omp_get_max_threads (); // upper bound may be distinct to whatever number user has chosen
#else
  n_threads = 1; // compiled without openMP support (e.g. --disable-openmp)
  biomcmc_warning ("Program compiled without multithread support");
#endif

  biomcmc_get_time (time0); 
  biomcmc_get_time (time1); 
  arg_parameters params = get_parameters_from_argv (argc, argv);
  if (params.dist->ival[0] < 0) params.dist->ival[0] = 0; 

  /* 1. read reference sequence (if missing, use database just to get length) */
  if (params.ref->count) {
    rfas = new_readfasta (params.ref->filename[0]);
    refseq = rfas->seq;
  }
  else {
    rfas = new_readfasta (params.fasta->filename[0]);
    refseq = NULL;
  }
  readfasta_next (rfas);
  
  if (params.trim->ival[0] > 0) trim = (size_t) params.trim->ival[0];
  if (trim > rfas->seqlength / 2.1) trim = rfas->seqlength / 2.1; // if we trim more than 1/2 the genome there's nothing left
  if (params.dist->ival[0] > (int)(rfas->seqlength) / 10) params.dist->ival[0] = rfas->seqlength / 10;

  /* 1.1 one cluster per thread, with char pointers etc*/
  clust = (cluster_t*) biomcmc_malloc (n_threads * sizeof (cluster_t));
  for (c = 0; c < n_threads; c++) clust[c] = new_cluster (refseq, rfas->seqlength, params.dist->ival[0], trim);
  seq_vec    = (char**) biomcmc_malloc (n_threads *sizeof (char*)); // only pointer
  name_vec   = (char**) biomcmc_malloc (n_threads *sizeof (char*)); // only pointer
  nchars_vec = (size_t*) biomcmc_malloc (n_threads *sizeof (size_t)); // actual number
  for (c = 0; c < n_threads; c++) { 
    seq_vec[c] = name_vec[c] = NULL; 
    nchars_vec[c] = 0;
  }

 
  del_readfasta (rfas);
  
  /* 2. read alignment files (can be several) */
  for (j = 0; j < params.fasta->count; j++) {
    rfas = new_readfasta (params.fasta->filename[j]);
    end_of_file = false;

    while (!end_of_file) {
#pragma omp single
      for (c = 0; c < n_threads; c++) { // reads n_threads sequences
        if (readfasta_next (rfas) >= 0) {
          seq_vec[c] = rfas->seq; rfas->seq = NULL;
          name_vec[c] = rfas->name; rfas->name = NULL;
          nchars_vec[c] = rfas->seqlength;
        }
        else {
          seq_vec[c] = name_vec[c] = NULL;
          end_of_file = true;
        }
      } // single thread for() loop

#pragma omp parallel for shared(c, count, clust, seq_vec, name_vec, nchars_vec, time1)
      for (c = 0; c < n_threads; c++) if (seq_vec[c]) {
        check_seq_against_cluster (clust[c], &(seq_vec[c]), &(name_vec[c]), nchars_vec[c]);
        if (!(++count % print_interval)) {
          fprintf (stderr, "Pool %3d has %d clusters; %d sequences analysed in total; last %d sequences took %4.4lf secs\n", 
                   c, clust[c]->n_fs, count, print_interval, biomcmc_update_elapsed_time (time1)); 
          fflush(stderr);
        }
      }
    } // while not end of file

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
 
  for (c = 0; c < n_threads; c++) qsort (clust[c]->fs, clust[c]->n_fs, sizeof (fastaseq_t), compare_fastaseq);

  // FIXME this assumes one cluster struct
  strcpy (outfilename + outlength, ".aln.xz");
  save_cluster_to_xz_file (clust[0], outfilename);

  strcpy (outfilename + outlength, ".csv.xz");
  save_neighbours_to_xz_file (clust[0], outfilename);
    
  fprintf (stderr, "Finished sorting clusters and saving files in %lf secs\n", biomcmc_update_elapsed_time (time0)); fflush(stderr);

  /* everybody is free to feel good */
  del_arg_parameters (params);
  for (c = 0; c < n_threads; c++) del_cluster (clust[c]); 
  if (clust) free (clust);
  if (seq_vec) free (seq_vec);
  if (name_vec) free (name_vec);
  if (nchars_vec) free (nchars_vec);
  if (outfilename) free (outfilename);
  return EXIT_SUCCESS;
}
