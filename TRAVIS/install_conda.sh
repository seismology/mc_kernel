#!/bin/bash

TRAVIS_ROOT="$1"

# Install Miniconda
export MINICONDA=$TRAVIS_ROOT/miniconda
export PATH="$MINICONDA/bin:$PATH"
hash -r
# Install conda only if necessary
if [ ! -f "$MINICONDA/bin/conda" ]; then
  wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -f -p $MINICONDA
fi
conda update --yes conda
conda install --yes -c obspy nomkl obspy netcdf4 psutil python=3.5
conda install --yes -c SciTools cartopy
