command -v conda >/dev/null || {
    wget http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh -O miniconda.sh;
    bash miniconda.sh -b -f -p $MINICONDA;
    conda update --yes conda;
    conda remove --yes -n condaenv --all;
    conda create --yes -n condaenv python=3.5;
    conda install --yes -n condaenv pip;
    source activate condaenv;
    conda install --yes -c obspy obspy h5py python=3.5 netcdf4 
  }
