<img src="recipe/uvaia-text.png" height="100" alt="Uvaia">

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/quadram-institute-bioscience/tatajuba/blob/master/LICENSE)

## Introduction

Uvaia is a program for pairwise reference-based alignment, and subsequent search against an aligned database. 
The alignment uses the promising [WFA library](https://github.com/smarco/WFA) implemented by Santiago Marco-Sola, and
the database search is based on score distances from my 
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) library.
In the past it used the [kseq.h](https://github.com/lh3/seqtk) library, by Heng Li, for reading fasta files, but
currently it relies on general compression libraries available on 
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib).
In particular I am trying to make all functions work with XZ compressed files.

Uvaia has been developed to help with SARS-CoV-2 analysis, and is been used as a quick replacement for [civet](https://github.com/COG-UK/civet).


#### Etymology
[Uvaia (Eugenia_pyriformis)](https://en.wikipedia.org/wiki/Eugenia_pyriformis) (also know as uvaieira, uaieira, ubaia e uvalha) 
is a fruit tree typical of Brazil. Its name comes from the tupi *iwa'ya*, which means "sour fruit", and if you have good 
imagination, its [pronunciation] resembles WFA.

## Installation
### Conda
[![Anaconda-Server Badge](https://anaconda.org/bioconda/uvaia/badges/platforms.svg)](https://anaconda.org/bioconda/uvaia)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/uvaia/badges/latest_release_date.svg)](https://anaconda.org/bioconda/uvaia)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/uvaia/badges/downloads.svg)](https://anaconda.org/bioconda/uvaia)

After you install [miniconda](https://conda.io/en/latest/miniconda.html), simply run
```[bash]
conda install -c bioconda uvaia
```
<!---
The version available in conda is outdated, please install it from source during the next few days.
If you really need a conda package, [here you can find](https://github.com/quadram-institute-bioscience/uvaia/issues/1#issuecomment-1092033238) 
a more recent version, but still outdated. 
-->

The conda version may not be up-to-date. The code is under active development while we prepare a manuscript for it. 

### Compiling from source
To install it from source, you should download this repository with `git clone --recursive` to ensure it also downloads its submodules (see below
for a **tl;dr**).
If you forgot to do so, you can update it with
```
git submodule update --init --recursive
```

It will compile from the submodules [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) 
and [our modified WFA](https://github.com/leomrtns/WFA) before finally compiling `uvaia`.
Notice that executables for other software it relies on are **not** generated, only their libraries are used as dependencies.

This sofware uses `autotools`, which means you can install it with `configure` + `make`.
You may need to define where you want it installed with `configure --prefix=DIR`, which is where are your unix-like
`include/`, `lib/`, and `bin/` directories. My favourite is `~/local`. Another popular option if you are on a conda
environment is `--prefix=${CONDA_PREFIX}`

Here is an example of its installation, please modify accordingly. 

```bash
/home/simpson/$ git clone --recursive git@github.com:quadram-institute-bioscience/uvaia.git
/home/simpson/$ cd uvaia && ./autogen.sh
/home/simpson/$ mkdir build && cd build
/home/simpson/$ ../configure --prefix=${HOME}/local
/home/simpson/$ make; make install
/home/simpson/$ make check  # battery of unit and integration tests, not mandatory
```
Remember that the installation and `autogen.sh` in particular modify/add local files; therefore updating the repository
from github will complain about uncommited changes. You can run `git stash` (or reinstall from scratch) before `git pull`.
check the directory [recipe](recipe/) for having a better idea of how conda/docker install it. 

If the compilation is unsuccessful, you can check if all libraries and packages below are installed:
```bash
## packages necessary for autotools, o.w. it will complain when you run "autogen.sh": 
/home/simpson/$ apt-get install pkg-config autotools-dev autoconf automake libtool
/home/simpson/$ (cd tatajuba && ./autogen.sh)  ## the parentheses avoid entering the directory afterwards

## C libraries needed or suggested by uvaia : 
/home/simpson/$ apt-get install zlib1g-dev libomp-dev libbz2-dev check liblzma-dev
/home/simpson/$ (cd build && ../configure)  ## etc.  
```

The libraries rely on `pkg-config` to find their location, otherwise you'll see a cryptic error message about `possibly
undefined macro`. 
If you can only install `pkg-config` through conda then you may need to install the C libraries via conda as well.
Or checking and updating your [`$PKG_CONFIG_PATH` environment variable](https://people.freedesktop.org/~dbn/pkg-config-guide.html).
Thus to install all dependent packages through conda (or mamba, which is faster :wink:):

```
mamba install automake libtool pkg-config make libgcc-ng check zlib xz bzip2 libgomp
```

## Running

The two main programs are:

* *uvaialign*, to align your query sequences against a reference genome. Output goes to `stdout` (your screen).
* *uvaia*, to search for _aligned_ queries against an _aligned_ database. Both query and database fasta files 
can be aligned to same reference sequence with `mafft`, `viralMSA`, or `uvaialign`. 
Outputs a file with a table of reference sequences which are closest neighbours to each query. 
It also generates a fasta file which includes these closest reference sequences (and also other, similar sequences).

There are also two experimental programs *uvaiaclust* and *uvaiaball*, which are still under development. 
And we also include the original `uvaia` program, called *uvaia_legacy* and used before 2022, which however cannot cope
with the huge file sizes we have currently.

### uvaialign
```
Align query sequences against a reference
The complete syntax is:

 uvaialign  [-h|--help] [-v|--version] [-a|--ambiguity=<double>] -r|--reference=<ref.fa|ref.fa.gz> <seqs.fa|seqs.fa.gz>

  -h, --help                       print a longer help and exit
  -v, --version                    print version and exit
  -a, --ambiguity=<double>         maximum allowed ambiguity for sequence to be excluded (default=0.5)
  -r, --reference=<ref.fa|ref.fa.gz> reference sequence
  <seqs.fa|seqs.fa.gz>             sequences to align
```
The mandatory arguments are the reference sequence (Wuhan-Hu-1 `MN908947.3` usually) and your unaligned sequences in
gzip or uncompressed format (xz format not yet available).

### uvaia

The main program, for searching for closest neighbours. Help page (called with `uvaia --help`):

```
For every query sequence, finds closest neighbours in reference alignment.
Notice that this software is multithreaded (and its performance depends on it)

The complete syntax is:

 uvaia  [-h|--help] [-v|--version] [--acgt] [-k|--keep_resolved] [-x|--exclude_self] [-n|--nbest=<int>] [-t|--trim=<int>] [-A|--ref_ambiguity=<double>] [-a|--query_ambiguity=<double>] [-p|--pool=<int>] -r|--reference=<ref.fa(.gz,.xz)> [-r|--reference=<ref.fa(.gz,.xz)>]... <seqs.fa(.gz,.xz)> [-t|--nthreads=<int>] [-o|--output=<without suffix>]

  -h, --help                       print a longer help and exit
  -v, --version                    print version and exit
  --acgt                           considers only ACGT sites (i.e. unambiguous SNP differences) in query sequences (mismatch-based)
  -k, --keep_resolved              keep more resolved and exclude redundant query seqs (default is to keep all)
  -x, --exclude_self               Exclude reference sequences with same name as a query sequence
  -n, --nbest=<int>                number of best reference sequences per query to store (default=100)
  -t, --trim=<int>                 number of sites to trim from both ends (default=0, suggested for sarscov2=230)
  -A, --ref_ambiguity=<double>     maximum allowed ambiguity for REFERENCE sequence to be excluded (default=0.5)
  -a, --query_ambiguity=<double>   maximum allowed ambiguity for QUERY sequence to be excluded (default=0.5)
  -p, --pool=<int>                 Pool size, i.e. how many reference seqs are queued to be processed in parallel (larger than number of threads, defaults to 64 per thread)
  -r, --reference=<ref.fa(.gz,.xz)> aligned reference sequences (can be several files)
  <seqs.fa(.gz,.xz)>               aligned query sequences
  -t, --nthreads=<int>             suggested number of threads (default is to let system decide; I may not honour your suggestion btw)
  -o, --output=<without suffix>    prefix of xzipped output alignment and table with nearest neighbour sequences
```

This is a rewrite of the old uvaia (currently available as `uvaia_legacy`) which works with very big reference
alignments (which cannot fit in memory at once).
It uses priority queues
(called `min_heap` internally) to store, for each query sequence, only the best scoring references. It assumes that the
database of reference sequences is too large to fit into memory, and thus we read the (possibly compressed) reference
alignment in batches, keeping only their names in memory and dumping all sequences that _at some point_ belonged to the
set. Thus in the beginning the software saves a lot of reference sequences which may not be very close to the query
sequences, and later it updates the output alignment less often, with closer sequences only. 

It runs in parallel using all available processors (in the future the user should be able to modify it), and it relies
on a "compression" of the query sequences into variable sites and common variants. It can also remove redundant
sequences, i.e. those identical to or equivalent but with fewer unambiguous sites with another.

#### output table

By default uvaia will generate an alignment file named `nn_uvaia.aln.xz` and a csv table `nn_uvaia.csv.xz`.
The csv table will look like
```
query,reference,rank,ACGT_matches,text_matches,partial_matches,valid_pair_comparisons,ACGT_matches_unique,valid_ref_sites
England/NORW-3078E97/2021,England/NORW-3061C36/2021,1,14985,14985,14987,14988,14985,29843
England/NORW-3078E97/2021,England/NORW-302EA07/2021,2,14984,14984,14986,14988,14984,29875
England/NORW-3078E97/2021,England/NORW-3034A7D/2021,3,14984,14984,14986,14988,14984,29869
England/NORW-3078E97/2021,England/NORW-3034B26/2021,4,14984,14984,14986,14988,14984,29851
England/NORW-3078E97/2021,England/NORW-306AB26/2021,5,14983,14983,14985,14988,14983,29813
England/NORW-3078E97/2021,England/NORW-31425AC/2022,6,14982,14982,14985,14988,14982,29169
```

Where the first three columns are the query sequence name, the reference sequence name, and the rank of the reference
with regard to the query &mdash; i.e. how close they are, such that rank=1 means that the reference is the closest to the
query, rank=2 means that it is the second closest to the query etc. The other columns are similarity (match) measures
and other numbers used in ranking them.
The output is ordered by rank, for each query. 

Unlike other distance calculation software, it actually calculates the number of matches between sequences as well as
the number of valid pairwise comparisons.
A distance can be created from the difference between `valid_pair_comparisons` (or genome length) and one of the number
of matches. 
As mentioned above, the rank (of how close the reference is to the query sequence) is given by the number of matches
etc, in the order given in the table below from the 4th to the 9th column.

column | column name | description
-----  | ----------- | -----------
1| *query*     | query sequence name
2| *reference* | reference sequence name (close to query sequence from column 1)
3| *rank* | order of closeness between reference and query, as given by columns 4-9
4| *ACGT_matches*        | considering only ACGT 
5| *text_matches*        |  exact matches, thus M-M is a match but M-A is not
6| *partial_matches*     | M-A is considered a match since the partially ambiguous `M` equals {A,C}. However the fully ambiguous `N` is neglected
7| *valid_pair_comparisons* | the `effective` sequence length for the comparison, i.e. sum of sites without gaps or N in any of the two sequences
8| *ACGT_matches_unique* | <nobr> a 'consensus' between query seqs is created, and this is the number of matches present in the query but not in the consensus (in short, it prefers neighbours farther from the common ancestor of the queries, in case of ties)  </nobr>
9| *valid_ref_sites*     | <nobr> if everything else is the same, then sequences with less gaps and Ns are preferred (caveat is that some sequencing labs artificially impute states, in practice removing all gaps and Ns) </nobr>
 
In other words, the rank is given by the number of `ACGT_matches`, and ties are broken by the number of `text_matches`,
and so forth.
The 4th, 5th, and 7th columns above are the most useful for the final user. But you can simply look at their rank, as
described below.  

#### output alignment

The alignment file will include all reference sequences temporarily assigned as closest to the queries. Thus it may have
many more sequences than the strict set of K nearest neighbours. 
If you want to extract the reference sequences which would compose the set of closest neighbours, you'll need to do something like:

```
xzcat nn_uvaia.csv.xz | cut -d "," -f 2 | sort | uniq > closest_names.txt
seqkit grep -f names.txt nn_uvaia.aln.xz > closest_names.aln
```

If you are using a recent version of `seqkit` (https://github.com/shenwei356/seqkit) which can handle XZ files.
By using the rank column from hte csv file, you can chose a subset of "best" sequences to add. For instance, to get only
the closest neighbours to each query sequence, the first command above could be replaced by

```
# reference sequence with rank = 1 (i.e. closest)
xzcat nn_uvaia.csv.xz | grep -e ",1," | cut -d "," -f 2 | sort | uniq > closest_names.txt

# all reference sequences with rank below 21 (i.e. the 20 closest references) 
xzcat nn_uvaia.csv.xz | gawk -F "," '$3 < 21 {print $2}' | sort | uniq > closest_names.txt
```

(Notice that the old version of uvaia, now called `uvaia_legacy`, created the final alignment with only the closest
references; this is not possible anymore in one pass due to very large sequence files). 

#### uvaia reports matches, not distances

The reported number of matches may differ between programs or runs due to how the query sequences are compressed and indexed, 
however their relative ranks should be preserved. 
Uvaia reports the total number of matches, which is a measure of similarity. 
Other programs report distances, which is a measure of dissimilarity. 

For instance `valid_pair_comparison` - `partial_matches` generates similar distances as [snp-dists](https://github.com/tseemann/snp-dists).
This is because [snp-dists](https://github.com/tseemann/snp-dists) excludes every non-ACGT, even partially informative sites, and counts only ACGT.
Here, sites with a gap or `N` in one of the sequences are ignored, while partially informative sites (e.g. `M` or `R`)
are included. 

The example sequences below illustrate the difference between them. All 3 sequences have 7 valid sites, since we exclude `N` and `-`:
```
seq1  AAC GTT A--    7 valid sites: 7 x ACGT + 0 partial
seq2  AAC G-T AM-    7 valid sites: 6 x ACGT + 1 partial (M)
seq3  MNC GTT MC-    7 valid sites: 5 x ACGT + 2 partial (M) 

ACGT_matches (seq1,seq2) = 6   partial_matches (seq1,seq2) = 6   valid_pair_comparisons(seq1,seq2) = 6
ACGT_matches (seq1,seq3) = 4   partial_matches (seq1,seq3) = 6   valid_pair_comparisons(seq1,seq3) = 6 
ACGT_matches (seq2,seq3) = 3   partial_matches (seq2,seq3) = 6   valid_pair_comparisons(seq2,seq3) = 6
```
So, although all sequences have a distance of zero from each other (no mismatches), they are not identical, and these
differences are relevant phylogenetically.
This highlights how `uvaia` and other software (in our example, snp-dists) see the sequences differently, remembering that only columns
where neither sequence has a `-` are considered:
```
      uvaia          snp-dists
seq1  AAC GTT A--    AAC GTT A--
seq2  AAC G-T AM-    AAC G-T A--
seq3  M-C GTT MC-    --C GTT -C-
```
I hope this example helps seeing why the difference between `valid_pair_comparison` and `partial_matches` (from uvaia)
is usually close to the number of ACGT mismatches (the default distances produced by `snp-dists`).
They do disagree if for instance *seq3* were instead `KNC GTT KC-`, since the snp-dists distance would be zero (`K`
is neglected, as is `M`), but uvaia knows that `K={G,T}` which is incompatible (and thus a mismatch) with `A` or `M={A,C}`.

#### matches are better for phylogenetic analysis

Usually uvaia should be used to find neighbouring sequences which will be used in downstream phylogenetic analysis.
Thus it keeps track of these several measures of similarity: in the limit, two very poor sequences with many gaps
and very few columns in common (e.g. `------AAA` and `CCC------`) have no mismatches!
The program periodically reports partial results, and it is not uncommon for the `Highest number of ACGT mismatches` to
increase as the run progresses. This apparent contradiction is in line with our objective of tracking **matches** instead of **mismatches**,
since a sequence pair with few mismatches may be due to too many ambiguous calls, and lower quality overall. And thus as
the program runs, it finds other sequences with many more _valid_ comparisons, at the cost of a few extra mismatches.
Furthermore the reported "highest number of mismatches" is the highest over query samples and neighbour sets, and 
may be due to one or a few query sequences.

However, **if you do want to reproduce other software behaviour and consider only As,Cs,Gs, and Ts**, there is an option
conveniently invoked as `--acgt`. With this option, it still tracks the matches but considering partially informative
sites (e.g. `M` or `K`) together with gaps and `N`s. The output table will be a bit different, with two extra columns
describing the distance (mismatches) from the reference to the consensus columns of the query, and the distance using
the "unique" (polymorphic) columns. If you are confused, don't worry and just use their sum as the "SNP distance". 
Or the difference between the number of valid comparisons and the number of matches, as usual.
If you use the option `--acgt`, then the generate files will be `nn_uvaia_acgt.aln.xz` and `nn_uvaia_acgt.csv.xz`.

### Other programs
Some description about the extra software included are provided in the [README_extra.md](README_extra.md) file, but I do
not advise you to use them.

## License 
SPDX-License-Identifier: GPL-3.0-or-later


Copyright (C) 2020-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

Tauari uses the library [WFA](https://github.com/smarco/WFA) by Santiago Marco-Sola, distributed under the MIT license.
It also uses [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) by Leonardo de Oliveira Martins
under a GPL3.0 license. 
The bearded goat from Uvaia's logo belongs to [Thomas Pennant's Allgemeine Uebersicht der vierfüssigen Thiere](https://www.flickr.com/photos/biodivlibrary/albums/72157715114535503) (public domain). 

![Anurag's github stats](https://github-readme-stats.vercel.app/api?username=leomrtns&count_private=true&show_icons=true&theme=calm)

[pronunciation]: https://www.google.com/search?sxsrf=ALeKk02RdD-T7ABhToHthWhPVrVOjEjrGQ:1600068291002&q=uvaia+pronounce&sa=X&ved=2ahUKEwi9u-fwjujrAhUqURUIHbzbAY0Q7xYoAHoECFoQKA&biw=1920&bih=1008
