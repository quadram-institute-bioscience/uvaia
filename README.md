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
It also uses the [kseq.h](https://github.com/lh3/seqtk) library, by Heng Li, for reading fasta files.

Uvaia has been developed to help with  SARS-CoV-2 analysis, and is been used as a quick replacement for [civet](https://github.com/COG-UK/civet).


#### Etymology
[Uvaia (Eugenia_pyriformis)](https://en.wikipedia.org/wiki/Eugenia_pyriformis) (also know as uvaieira, uaieira, ubaia e uvalha) 
is a fruit tree typical of Brazil. Its name comes from the tupi *iwa'ya*, which means "sour fruit", and if you have good 
imagination, its [pronunciation] resembles WFA.

## Installation
### Conda
[![Anaconda-Server Badge](https://anaconda.org/bioconda/tatajuba/badges/platforms.svg)](https://anaconda.org/bioconda/uvaia)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/tatajuba/badges/latest_release_date.svg)](https://anaconda.org/bioconda/uvaia)

After you install [miniconda](https://conda.io/en/latest/miniconda.html), simply run
```[bash]
conda install -c bioconda uvaia
```

The conda version may not be up-to-date, currently it works with gzip but not LZMA (`.xz`) files.

### Compiling from source
If you want to install it from source (e.g. to use the `.xz` format), you should download this repository with `git clone --recursive` to ensure it also downloads its submodules (see below
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

```[bash]
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
```[bash]
## packages necessary for autotools (o.w. it will complain when you run "autogen.sh" : 
/home/simpson/$ apt-get install pkg-config autotools-dev autoconf automake libtool
/home/simpson/$ (cd tatajuba && autoreconf)  ## the parentheses avoid entering the directory afterwards
## C libraries needed or suggested by uvaia : 
/home/simpson/$ apt-get install zlib1g-dev libomp-dev libbz2-dev check liblzma-dev
/home/simpson/$ (cd build && ../configure)  ## etc.  
```

## Running

The two main programs are:

* *uvaialign*, to align your query sequences against a reference genome. Output goes to `stdout` (your screen).
* *uvaia*, to search for _aligned_ queries against an _aligned_ database. Both query and database fasta files 
can be aligned to same reference sequence with `mafft`, `viralMSA`, or `uvaialign`. Output is a table of reference
sequences which are closest neighbours to each query. It can also generate a fasta file with these reference sequences. 

There is also an experimental program *uvaiaclust* which tries to remove redundant sequences from an alignment. 

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

Currently this program is a bit slow (10k sequences in one hour?).

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
