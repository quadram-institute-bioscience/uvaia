#!/bin/sh


if [ ! -e "./biomcmc-lib" ]; then  # -e is file (-f), directory (-d), link (-h) etc.
  ln -s submodules/biomcmc-lib ./biomcmc-lib

export AUTOMAKE="automake --foreign -a"
autoreconf -f -i
