#include <biomcmc.h>

#include "fastaseq.h"

typedef struct cqueue_struct* cqueue_t;

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_int  *dist;
  struct arg_int  *trim;
  struct arg_int  *pool;
  struct arg_int  *snps;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_file *out;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

struct cqueue_struct
{
  int n_clust;
  size_t *nchars;
  char **seq, **name;
  cluster_t *clust;
};

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

cqueue_t new_cqueue (int n_clust, char *refseq, size_t nsites, int dist, int trim, int snps);
char* get_outfile_prefix (const char *prefix, size_t *length);
char* read_reference_sequence (readfasta_t *rfas, const char *filename, int nseqs);
cqueue_t new_cqueue (int n_clust, char *refseq, size_t nsites, int dist, int trim, int snps);
void del_cqueue (cqueue_t cq);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .dist    = arg_intn("d","distance", NULL, 0, 1, "seqs with this SNP differences or less will be merged (default=1)"),
    .trim    = arg_int0(NULL,"trim", NULL, "number of sites to trim from both ends (default=0, suggested for sarscov2=230)"),
    .pool    = arg_int0("p","pool", NULL, "Pool size, i.e. number of clustering queues (should be larger than avail threads)"),
    .snps    = arg_int0("s","snps", NULL, "how many SNPs w.r.t. reference it keeps track (default=1, should be small number)"),
    .ref     = arg_filen("r","reference", "<ref.fa(.gz,.xz)>", 0, 1, "reference sequence (medoids are furthest from it)"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa(.gz,.xz)>", 1, 1024, "alignments to merge"),
    .out     = arg_file0("o","output", "<without suffix>", "prefix of xzipped output alignment and cluster table files"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.dist, params.trim, params.pool, params.snps, params.ref, params.fasta, params.out, params.end};
  params.argtable = argtable; 
  params.dist->ival[0] = 1;
  params.trim->ival[0] = 0;
  params.snps->ival[0] = 1;
#ifdef _OPENMP
  params.pool->ival[0] = 4 * omp_get_max_threads (); // default is to have quite a few
#else
  params.pool->ival[0] = 4;
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
  if (params.dist)  free (params.dist);
  if (params.trim) free (params.trim);
  if (params.pool) free (params.pool);
  if (params.snps) free (params.snps);
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
    printf ("Otherwise a medoid is the most resolved (less Ns, with no SNP contradicting previous medoid) in cluster\n.");
    printf ("A pool of independent clustering queues is created, such that each sequence is compared to only one of them at first.\n\n");
  }

  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, j, c, count = 0, n_clust = 256, print_interval = 5000, save_interval = 10000;
  bool end_of_file = false;
  int64_t time0[2], time1[2];
  double elapsed = 0.;
  size_t outlength = 0;
  char *outfilename = NULL, *refseq = NULL;
  readfasta_t rfas;
  cqueue_t cq;

  biomcmc_get_time (time0); 
  biomcmc_get_time (time1); 
  arg_parameters params = get_parameters_from_argv (argc, argv);
  if (params.dist->ival[0] < 0) params.dist->ival[0] = 0; 
  if (params.snps->ival[0] < 0) params.snps->ival[0] = 0;

  fprintf (stderr, "Experimental program: %s package: %s\n", basename(argv[0]), PACKAGE_STRING);
  
#ifdef _OPENMP
  n_clust = omp_get_max_threads (); // upper bound may be distinct to whatever number user has chosen
#else
  n_clust = 1; // compiled without openMP support (e.g. --disable-openmp)
  biomcmc_warning ("Program compiled without multithread support");
#endif
  if (params.pool->ival[0] >= n_clust) n_clust = params.pool->ival[0];
  fprintf (stderr, "Creating a pool of %d cluster queues; maximum distance is %d, and %d SNP locations are kept\n", 
           n_clust,params.dist->ival[0], params.snps->ival[0]); 
  fflush(stderr);

  
  if (params.out->count) outfilename = get_outfile_prefix (params.out->filename[0], &outlength);
  else                   outfilename = get_outfile_prefix ("cluster_uvaia", &outlength);

  /* 1. read reference sequence. If missing, use database to create a dummy sequence (composed of ACGT) */
  if (params.ref->count) refseq = read_reference_sequence (&rfas, params.ref->filename[0], 1);
  else                   refseq = read_reference_sequence (&rfas, params.fasta->filename[0], 1024);
  

  /* 1.1 one cluster per thread, with char pointers etc*/
  cq = new_cqueue (n_clust, refseq, rfas->seqlength, params.dist->ival[0], params.trim->ival[0], params.snps->ival[0]);
 
  if (refseq) free (refseq);
  del_readfasta (rfas);
  
  /* 2. read alignment files (can be several) and fill pool of cluster queues */
  count = 0;
  for (j = 0; j < params.fasta->count; j++) {
    rfas = new_readfasta (params.fasta->filename[j]);
    end_of_file = false;

    while (!end_of_file) {
#pragma omp single
      for (c = 0; c < n_clust; c++) { // reads n_clust sequences
        if (readfasta_next (rfas) >= 0) {
          count++;
          cq->seq[c] = rfas->seq; rfas->seq = NULL;
          cq->name[c] = rfas->name; rfas->name = NULL;
          cq->nchars[c] = rfas->seqlength;
        }
        else {
          cq->seq[c] = cq->name[c] = NULL;
          end_of_file = true;
        }
      } // single thread for() loop

#pragma omp parallel for shared(c, count, cq)
      for (c = 0; c < cq->n_clust; c++) if (cq->seq[c]) {
        check_seq_against_cluster (cq->clust[c], &(cq->seq[c]), &(cq->name[c]), cq->nchars[c]);
        if (!(cq->clust[c]->n_fs % print_interval) && !(count%4)) { // cont%4 to avoid consecutive prints
          fprintf (stderr, "Queue %3d / %d has %d clusters; %d sequences analysed in total;\n", c, cq->n_clust, cq->clust[c]->n_fs, count); 
          fflush(stderr);
        }
      }

      if ((count >= print_interval) && ((count % print_interval) < n_clust)) {
        elapsed = biomcmc_update_elapsed_time (time1); 
        fprintf (stderr, "%d sequences analysed in total; last %d sequences took %.3lf secs\n", count, print_interval, elapsed); 
        fflush(stderr);
      }
      if ((count >= save_interval) && ((count % save_interval) < n_clust) && (elapsed > 30)) {
        fprintf (stderr, "Saving partial clustering info from %d sequences to file %s\n", count, outfilename); fflush(stderr);
        save_neighbours_to_xz_file (cq->clust, cq->n_clust, outfilename);
      }
    } // while not end of file

    del_readfasta (rfas);
    fprintf (stderr, "Finished reading file %s in %.3lf secs; Commulative %d sequences read\n", params.fasta->filename[j], biomcmc_update_elapsed_time (time0), count); fflush(stderr);
  }  // for fasta file

  /* 2.1 create vector of SNP locations wrt reference */
  generate_idx_from_cluster_list (cq->clust, cq->n_clust, 0);

//#pragma omp parallel for shared(c, clust, j) private (count)
//  for (c = 0; c < n_clust; c++) {
//   count = compact_cluster (clust[c]);  // does not improve at all (zero coalescences)
//    fprintf (stderr, "%d clusters coalesced within queue %d\n", count, c); fflush(stderr);
//  }

  /* 3. merge clusters */

  elapsed = biomcmc_update_elapsed_time (time1); // time1 is more finegrained than time0 
  for (c = cq->n_clust; c > 1; c = (c/2 + c%2)) { // reduce (outside parallel loop)
#pragma omp parallel for shared(cq, c, j) private(count,i)
    for (j = 0; j < c/2; j++) {
      i = j + c/2 + c%2;
      fprintf (stderr, "Merging queues %d and %d out of %d (%d and %d nonredundant sequences)\n",j,i,c,cq->clust[j]->n_fs, cq->clust[i]->n_fs); fflush(stderr);
      count = merge_clusters (cq->clust[j], cq->clust[i]); 
      fprintf (stderr, "%d clusters coalesced between queues %d and %d\n", count, j, i); fflush(stderr);
    }

    fprintf (stderr, "%.3lf secs elapsed. Saving partial clustering info from %d sequences to file %s\n", biomcmc_update_elapsed_time (time1), count, outfilename); fflush(stderr);
    save_neighbours_to_xz_file (cq->clust, cq->n_clust, outfilename);
  }

  //for (c = 0; c < n_clust; c++) qsort (clust[c]->fs, clust[c]->n_fs, sizeof (fastaseq_t), compare_fastaseq);
  qsort (cq->clust[0]->fs, cq->clust[0]->n_fs, sizeof (fastaseq_t), compare_fastaseq);

  save_neighbours_to_xz_file (cq->clust, 1, outfilename);
  strcpy (outfilename + outlength, ".aln.xz");
  save_cluster_to_xz_file (cq->clust, 1, outfilename);

  fprintf (stderr, "Finished sorting clusters and saving files in %lf secs\n", biomcmc_update_elapsed_time (time0)); fflush(stderr);

  /* everybody is free to feel good */
  del_arg_parameters (params);
  del_cqueue (cq);
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
  strcpy (filename + *length, ".csv.xz"); // CSV is saved several times
  return filename;
}
  
char*
read_reference_sequence (readfasta_t *rfas, const char *filename, int nseqs)
{
  char *refseq = NULL;
  int j = -1, i, count = 0xff; // count is number of non-ACGT sites
  *rfas = new_readfasta (filename);
  
  fprintf (stderr, "Generating a reference from up to %d sequences in %s\n", nseqs, filename);
  for (i = 0; (readfasta_next (*rfas) > 0) && (count) && (i < nseqs); i++) {
    if (j < 0) j = (int) (*rfas)->seqlength;
    else if (j != (int) (*rfas)->seqlength) biomcmc_error ("Unaligned sequences: first seq has %d sites but %s has %lu sites\n", j, (*rfas)->name, (*rfas)->seqlength);
    count = accumulate_reference_sequence (&refseq, (*rfas)->seq, (*rfas)->seqlength); // merging of ACGT sites
  }
  count = replace_Ns_from_reference (refseq, (*rfas)->seqlength);
  if (count) fprintf (stderr, "Reference still had %d Ns or indels which were replaced arbitrarily (it only makes program a bit slower).", count);
  return refseq;

}

cqueue_t
new_cqueue (int n_clust, char *refseq, size_t nsites, int dist, int trim, int snps)
{
  cqueue_t cq = (cqueue_t) biomcmc_malloc (sizeof (struct cqueue_struct));
  int c;
  cq->n_clust = n_clust;
  cq->clust = (cluster_t*) biomcmc_malloc (n_clust * sizeof (cluster_t));

  if (trim < 0) trim = 0;
  if (trim > nsites / 2.1) trim = nsites / 2.1; // if we trim more than 1/2 the genome there's nothing left
  if (dist > (int)(nsites) / 10) dist = (int) nsites / 10;

  for (c = 0; c < n_clust; c++) cq->clust[c] = new_cluster (refseq, nsites, dist, (size_t) trim, snps);
  cq->seq    = (char**) biomcmc_malloc (n_clust *sizeof (char*)); // only pointer
  cq->name   = (char**) biomcmc_malloc (n_clust *sizeof (char*)); // only pointer
  cq->nchars = (size_t*) biomcmc_malloc (n_clust *sizeof (size_t)); // actual number
  for (c = 0; c < n_clust; c++) { 
    cq->seq[c] = cq->name[c] = NULL; 
    cq->nchars[c] = 0;
  }
  return cq;
}

void
del_cqueue (cqueue_t cq)
{
  if (!cq) return;
  int c;
  for (c = cq->n_clust -1; c >= 0; c--) del_cluster (cq->clust[c]); 
  if (cq->clust) free (cq->clust);
  if (cq->seq) free (cq->seq);
  if (cq->name) free (cq->name);
  if (cq->nchars) free (cq->nchars);
  free (cq);
  return;
}
