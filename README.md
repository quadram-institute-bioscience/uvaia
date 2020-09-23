<img src="recipe/uvaia-text.png" height="100" alt="Uvaia">

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Introduction

Uvaia is an experimental program for pairwise reference-based alignment, and subsequent search against an aligned database. 
The alignment uses the promising [WFA library](https://github.com/smarco/WFA) implemented by Santiago Marco-Sola, and
the database search is based on score distances from my 
[biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) library.
It also uses the [kseq.h](https://github.com/lh3/seqtk) library, by Heng Li, for reading fasta files.

Uvaia has been developed to help with  SARS-CoV-2 analysis.

(this is still a weekends project so progress may be slow.)

#### Etymology
[Uvaia (Eugenia_pyriformis)](https://en.wikipedia.org/wiki/Eugenia_pyriformis) (also know as uvaieira, uaieira, ubaia e uvalha) 
is a fruit tree typical of Brazil. Its name comes from the tupi *iwa'ya*, which means "sour fruit", and if you have good 
imagination, its [pronunciation] resembles WFA.

## Installation
You should download this repository with `git clone --recursive` to ensure it also downloads its submodules (see below
for a **tl;dr**).
If you forgot to do so, you can update it with
```
git submodule update --init --recursive
```

It will compile from the submodules `biomcmc-lib` and `WFA` before finally compiling `uvaia`.
Notice that executables for other software it relies on are **not** generated, only their libraries are used as dependencies.

This sofware uses `autotools`, which means you can install it with `configure` + `make`.
You may need to define where you want it installed with `configure --prefix=DIR`, which is where are your unix-like
`include/`, `lib/`, and `bin/` directories. My favourite is `~/local`. Another popular option if you are on a conda
environment is `--prefix=${CONDA_PREFIX}`

Here is an example of its installation, please modify accordingly. 

```[bash]
/home/simpson/$ git clone --recursive git@github.com:leomrtns/uvaia.git
/home/simpson/$ mkdir build && cd build
/home/simpson/$ ../uvaia/configure --prefix=${HOME}/local
/home/simpson/$ make; make install
/home/simpson/$ make check  # battery of unit and integration tests, not mandatory
```

## Running

This is a very beta project so usage may change overnight... in any case, two programs are emerging:

* *uvaialign*, to align your query sequences against a reference genome. Output goes to `stdout` (your screen).
* *uvaia*, to search for _aligned_ queries against an _aligned_ database. Both query and database fasta files 
can be aligned to same reference sequence with `mafft`, `viralMSA`, or `uvaialign`. Output is a table of reference
sequences which are closest neighbours to each query. It can also generate a fasta file with these reference sequences. 

## License 
SPDX-License-Identifier: GPL-3.0-or-later

Copyright (C) 2020-today  [Leonardo de Oliveira Martins](https://github.com/leomrtns)

This is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version (http://www.gnu.org/copyleft/gpl.html).

Tauari uses the library [WFA](https://github.com/smarco/WFA) by Santiago Marco-Sola, distributed under the MIT license.
It also uses [biomcmc-lib](https://github.com/quadram-institute-bioscience/biomcmc-lib) by Leonardo de Oliveira Martins
under a GPL3.0 license. 

![Anurag's github stats](https://github-readme-stats.vercel.app/api?username=leomrtns&count_private=true&show_icons=true&theme=calm)

[pronunciation]: https://www.google.com/search?sxsrf=ALeKk02RdD-T7ABhToHthWhPVrVOjEjrGQ:1600068291002&q=uvaia+pronounce&sa=X&ved=2ahUKEwi9u-fwjujrAhUqURUIHbzbAY0Q7xYoAHoECFoQKA&biw=1920&bih=1008
