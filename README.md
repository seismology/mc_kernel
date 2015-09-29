MC Kernel
======

Calculate seismic sensitivity kernels using wavefields created with Axisem.

Authors:
Simon Stähler

Martin van Driel

Ludwig Auer

Kasra Hosseini

Tarje Nissen-Meyer

## Download
To download the current developer version, use git clone
```
git clone https://github.com/sstaehler/kerner
```

## Prequisites
MC Kernel needs a recent Fortran compiler (gfortran 4.8 or later), an MPI installation, the NetCDF library for accessing the wavefield files, the FFTW library for Fourier transforms and LAPACK, because we were too lazy to write a matrix inversion routine ourselves. The installation of these libraries can be done by hand on desktop machines or using modules on HPC environments
#### Ubuntu/Debian Linux
Since Ubuntu 14.04LTS, the system libraries can be used:
```bash
sudo apt-get install gfortran libnetcdff5 libfftw3-dev libblas3 openmpi-bin gfortran
```

#### MacOS X
We recommend using Homebrew to download and compile the necessary libraries. Be careful to add the *--with-fortran* argument while installing any library. Otherwise, only C libraries are installed. If you did install the libraries before without the *--with-fortran* argument, it may be necessary to remove and reinstall them.

#### HPC environments
Be sure to load modules for Fortran, NetCDF, FFTW and LAPACK.
On SuperMUC, the necessary commands are
```
module load fortran
module load netcdf/4.2
module load fftw
module load mkl
```

## Installation
Change into the download directory and copy the included template files into the main directory
```bash
cd kerner
./copy_templates.sh
```
The file *make_mc_kerner.macros* allows you to modify the compiler name and compiler flags according to your system. To compile, use make
```bash
make 
```
and 
```bash
make test 
```
to run a set of tests.

## FAQs
### Anaconda interference 
The python package manager Anaconda installs his own version of gfortran, which interferes with the system version and the MPI libraries. 
A solution is to replace the compiler name *FC* and the *MPIRUN* variable in *make_mc_kerner.macros* with absolute paths
```bash
FC = /usr/bin/mpif90
```
