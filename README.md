<img src="recipe/.png" height="100">

__Leonardo de Oliveira Martins<sup>1</sup>__
<br>
<sub>1. Quadram Institute Bioscience, Norwich Research Park, NR4 7UQ, UK</sub>

## Introduction

mussum dorma

## Installation
You should download this repository with `git clone --recursive` to ensure it also downloads its submodules.
It will compile from the directories `biomcmc-lib` and `WFA` before finally compiling `uvaia`.
Notice that executables for libraries above are **not** generated, only their libraries are used as dependencies.

This sofware uses `autotools`, so you can install it with `configure` + `make`.
You may need to define where you want it installed with `configure --prefix=DIR`, which is where are your unix-like
`include/`, `lib/`, and `bin/` directories. My favourite is `~/local`.

Here is an example of its installation, please modify accordingly. 

```[bash]
/home/simpson/$ git clone --recursive git@github.com:leomrtns/uvaia.git
/home/simpson/$ mkdir build && cd build
/home/simpson/$ ../uvaia/configure --prefix=${HOME}/local
/home/simpson/$ make; make install
/home/simpson/$ make check  # battery of unit and integration tests, not mandatory
```

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

