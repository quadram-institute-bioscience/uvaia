#include <biomcmc.h>

#include "fastaseq.h"

typedef struct queue_struct* queue_t;

typedef struct
{
  struct arg_lit  *help;
  struct arg_lit  *version;
  struct arg_int  *dist;
  struct arg_int  *trim;
  struct arg_int  *pool;
  struct arg_file *ref;
  struct arg_file *fasta;
  struct arg_file *out;
  struct arg_end  *end;
  void **argtable;
} arg_parameters;

struct queue_struct
{
  int n_seqs, *mindist;
  char **seq, **name; // should have a mindist for query seqs, for each element of queue to avoid race conditions
};

arg_parameters get_parameters_from_argv (int argc, char **argv);
void del_arg_parameters (arg_parameters params);
void print_usage (arg_parameters params, char *progname);

char* get_outfile_prefix (const char *prefix, size_t *length);

queue_t new_queue (int n_seqs, alignment query, int dist, int trim);
void del_queue (queue_t cq);

arg_parameters
get_parameters_from_argv (int argc, char **argv)
{
  arg_parameters params = {
    .help    = arg_litn("h","help",0, 1, "print a longer help and exit"),
    .version = arg_litn("v","version",0, 1, "print version and exit"),
    .dist    = arg_intn("d","distance", NULL, 0, 1, "ball radius, i.e. refs within this distance to any query seq are kept"),
    .trim    = arg_int0("t","trim", NULL, "number of sites to trim from both ends (default=0, suggested for sarscov2=230)"),
    .pool    = arg_int0("p","pool", NULL, "Pool size, i.e. how many reference seqs are queued to be processed in parallel (can be much larger than avail threads, defaults to 100)"),
    .ref     = arg_filen("r","reference", "<ref.fa(.gz,.xz)>", 1, 1024, "aligned reference sequences (can be several files)"),
    .fasta   = arg_filen(NULL, NULL, "<seqs.fa(.gz,.xz)>", 1, 1, "aligned query sequences"),
    .out     = arg_file0("o","output", "<without suffix>", "prefix of xzipped output alignment with subset of ref sequences"),
    .end     = arg_end(10) // max number of errors it can store (o.w. shows "too many errors")
  };
  void* argtable[] = {params.help, params.version, params.dist, params.trim, params.pool, params.snps, params.ref, params.fasta, params.out, params.end};
  params.argtable = argtable; 
  params.dist->ival[0] = 5;
  params.trim->ival[0] = 0;
#ifdef _OPENMP
  params.pool->ival[0] = 100 * omp_get_max_threads (); // default is to have quite a few
#else
  params.pool->ival[0] = 100;
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
  printf ("finds all reference neighbours to query sequences within a ball radius (i.e. close to any query sequence)\n");
  printf ("The complete syntax is:\n\n %s ", basename(progname));
  arg_print_syntaxv (stdout, params.argtable, "\n\n");
  arg_print_glossary(stdout, params.argtable,"  %-32s %s\n");
  if (params.help->count) {
    printf ("experimental algorithm; if radius smaller than closest neighour to a sequence, then this sequence will have no neighbours!\n\n");
  }

  del_arg_parameters (params);
  if (params.end->count && (!params.help->count)) exit (EXIT_FAILURE);
  exit (EXIT_SUCCESS);
}

int
main (int argc, char **argv)
{
  int i, j, c, count = 0, n_clust = 256, print_interval = 50000;
  bool end_of_file = false;
  int64_t time0[2], time1[2];
  double elapsed = 0.;
  size_t outlength = 0;
  char *outfilename = NULL, *refseq = NULL;
  alignment query;  // all query sequences in memory 
  readfasta_t rfas; // single fasta sequence (for queue of ref seqs)
  queue_t cq;       // queue of ref sequences

  biomcmc_get_time (time0); 
  arg_parameters params = get_parameters_from_argv (argc, argv);
  if (params.dist->ival[0] < 0) params.dist->ival[0] = 0; 

#ifdef _OPENMP
  n_clust = omp_get_max_threads (); // upper bound may be distinct to whatever number user has chosen
#else
  n_clust = 1; // compiled without openMP support (e.g. --disable-openmp)
  biomcmc_warning ("Program compiled without multithread support");
#endif
  if (params.pool->ival[0] >= n_clust) n_clust = params.pool->ival[0];
  fprintf (stderr, "Creating a queue of %d sequences; radius distance is %d\n", n_clust,params.dist->ival[0]); 
  fflush(stderr);

  
  if (params.out->count) outfilename = get_outfile_prefix (params.out->filename[0], &outlength); // outfilename already has "aln.xz" suffix
  else                   outfilename = get_outfile_prefix ("ball_uvaia", &outlength);

  /* 1. read query sequences into char_vectors */
  query = read_fasta_alignment_from_file ((char*)params.fasta->filename[0], 0xf); // 0xf -> neither true or false, but gets only bare alignment info

  if (params.trim->ival[0] <= 0) trim = 0; // lower bound is zero (default value)
  else trim = (size_t) params.trim->ival[0];
  if (trim > query->nchar / 2.1) trim = query->nchar / 2.1; // if we trim more than 1/2 the genome there's nothing left

  fprintf (stderr, "Finished reading %d query references in %lf secs;\n", query->ntax, biomcmc_update_elapsed_time (time0)); fflush(stderr);
  biomcmc_get_time (time1); // starts counting here (since we'll have a partial time0 and a total time1 chronometer 

  /* 1.1 several ref sequences per thread (total n_clust read at once), with char pointers etc*/
  cq = new_queue (n_clust, query, params.dist->ival[0], params.trim->ival[0]);
 
  /* 2. read alignment files (can be several) and fill pool of cluster queues */
  count = 0;
  for (j = 0; j < params.ref->count; j++) {
    rfas = new_readfasta (params.ref->filename[j]);
    end_of_file = false;

    while (!end_of_file) {
#pragma omp single
      for (c = 0; c < n_clust; c++) { // reads n_clust sequences
        if (readfasta_next (rfas) >= 0) {
          count++;
          cq->seq[c] = rfas->seq; rfas->seq = NULL;
          cq->name[c] = rfas->name; rfas->name = NULL;
          if (rfas->seqlength != query->nchar) {
            biomcmc_warning ("Reference sequence '%s' has %d sites but query sequences have %d sites\n", cq->name[c], rfas->seqlength, query->nchar);
            del_queue (cq);
            del_alignment (query);
            del_readfasta (rfas);
            biomcmc_error ("all sequences must be aligned");
          }
        }
        else {
          cq->seq[c] = cq->name[c] = NULL;
          end_of_file = true;
        }
      } // single thread for() loop

#pragma omp parallel for shared(c, count, cq)
      for (c = 0; c < cq->n_clust; c++) if (cq->seq[c]) {
        cq->mindist[c] = 0xffff;
        seq_ball_against_alignment (&(cq->seq[c]), &(cq->name[c]), &(cq->mindist[c]), query);
      }

      if ((count >= print_interval) && ((count % print_interval) < n_clust)) {
        elapsed = biomcmc_update_elapsed_time (time1); 
        fprintf (stderr, "%d sequences analysed in total; last %d sequences took %.3lf secs\n", count, print_interval, elapsed); 
        fflush(stderr);
      }

#pragma omp single
      for (c = 0; c < n_clust; c++) { // write seqs to fasta FIXME 
        //save_cluster_to_xz_file (cq->clust, 1, outfilename);
      }
    } // while not end of file

    del_readfasta (rfas);
    fprintf (stderr, "Finished reading file %s in %.3lf secs; Cummulative %d sequences read\n", params.ref->filename[j], biomcmc_update_elapsed_time (time0), count); fflush(stderr);
  }  // for fasta file

  /* everybody is free to feel good */
  del_arg_parameters (params);
  del_queue (cq);
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
  strcpy (filename + *length, ".aln.xz");  
  return filename;
}
  
queue_t
new_queue (int n_seqs, alignment query, int dist, int trim)
{
  queue_t cq = (queue_t) biomcmc_malloc (sizeof (struct queue_struct));
  int c;
  cq->n_seqs = n_seqs;

  if (trim < 0) trim = 0;
  if (trim > nsites / 2.1) trim = nsites / 2.1; // if we trim more than 1/2 the genome there's nothing left
  if (dist > (int)(nsites) / 10) dist = (int) nsites / 10;

  cq->seq    = (char**) biomcmc_malloc (n_seqs * sizeof (char*)); // only pointer, actual seq is allocated by readfasta
  cq->name   = (char**) biomcmc_malloc (n_seqs * sizeof (char*)); // only pointer '' '' ''
  cq->mindist = (int*)  biomcmc_malloc (n_seqs * sizeof (int*));
  for (c = 0; c < n_seqs; c++) { cq->seq[c] = cq->name[c] = NULL; cq->mindist[c] = 0xffff; }
  return cq;
}

void
del_queue (queue_t cq)
{
  if (!cq) return;
  int c;
  if (cq->seq) free (cq->seq);
  if (cq->name) free (cq->name);
  free (cq);
  return;
}
