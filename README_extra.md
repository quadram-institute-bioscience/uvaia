# Extra programs (use at your own risk)

This are unfinished or superceded programs. I will probably do not take action for issues reported relating to them. 

### uvaia_legacy

This is the original version of `uvaia`, which is mentioned in papers before 2022. The results can be reproduced using
the new version, by the way, but I'll leave the original program and code here for reproducibility. 

```
Search query sequences against reference ones to describe closest ones
This program relies on aligned sequences (queries and references must have same size)
The complete syntax is:

 uvaia  [-h|--help] [-v|--version] [-n|--nbest=<int>] [-m|--nmax=<int>] [--trim=<int>] [-A|--ref_ambiguity=<double>] [-a|--query_ambiguity=<double>] -r|--reference=[ref.fa(.gz)] [-o|--output=[chosen_refs.fa.xz]] [-t|--nthreads=<int>] [query.fa(.gz,.xz)]

  -h, --help                       print a longer help and exit
  -v, --version                    print version and exit
  -n, --nbest=<int>                number of best reference sequences per query to show (default=8)
  -m, --nmax=<int>                 max number of best reference sequences when several optimal (default=2 x nbest)
  --trim=<int>                     number of sites to trim from both ends (default=0, suggested for sarscov2=230)
  -A, --ref_ambiguity=<double>     maximum allowed ambiguity for REFERENCE sequence to be excluded (default=0.5)
  -a, --query_ambiguity=<double>   maximum allowed ambiguity for QUERY sequence to be excluded (default=0.5)
  -r, --reference=[ref.fa(.xz,.gz,.bz)]    *aligned* reference sequences
  -o, --output=[chosen_refs.fa.xz] XZIPPED (LZMA) output reference sequences (default is to not save sequences)
  -t, --nthreads=<int>             suggested number of threads (default is to let system decide; I may not honour your suggestion btw)
  [query.fa(.gz,.xz)]              *aligned* sequences to search for neighbour references
```

Both the database of aligned sequences and the set of aligned query sequences can be compressed files (xz, bz, gz). 
This program uses a lot of memory since it stores the whole (uncompressed) database in memory.

### uvaiaclust (experimental)

Removal of redundant sequences based on a single-distance canopy clustering: each new sequence will be merged into the
first cluster s.t. its distance is smaller than the threshold. Each cluster is represented by its most resolved sequence
(fewer Ns). An aligment with the final medoids is returned, with a list of sequence names belonging to each cluster. 

It is not a proper clustering since identical sequences may belong to distinct medoids. 
Furthermore the pairwise distances within the cluster may exceed the threshold (since medoids are updated as new
elements are added). 

If a reference sequence is provided, it is used when updating the medoids, which should be the most resolved and
furthest from the reference. 
This is to penalise assemblies where a low resolution is masqueraded by imputing the reference base.

Currently this program is a bit slow (10k sequences in one hour?), and some improvements used in uvaia haven't been added
here.

### uvaiaball (experimental)

This is a ball radius-based neighbour search, but doesn't have some improvements added to `uvaia`
