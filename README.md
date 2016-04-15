# MC Kernel
[![DOI](https://zenodo.org/badge/21468/seismology/mc_kernel.svg)](https://zenodo.org/badge/latestdoi/21468/seismology/mc_kernel)
![](https://www.geophysik.uni-muenchen.de/~staehler/kerner/logo.png)

### Authors:
Simon Stähler, Martin van Driel, Ludwig Auer, Kasra Hosseini, Tarje Nissen-Meyer

### License: 
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, version 3 of the License.

### Acknowledgements:
If you use this program for your work, please cite:

Simon C. Stähler, Martin van Driel, Ludwig Auer, Kasra Hosseini, Karin Sigloch, and Tarje Nissen-Meyer: **MC Kernel: Broadband Waveform Sensitivity Kernels for Seismic Tomography**, *Geophysical Research Abstracts*, Vol. 18, EGU2016-7020-2, 2016

### Download
To download the current developer version, use git clone
```
git clone https://github.com/seismology/mc_kernel
```

## Prequisites
MC Kernel needs a recent Fortran compiler (gfortran 4.8 or later), an MPI installation, the NetCDF library for accessing the wavefield files, the FFTW library for Fourier transforms and LAPACK, because we were too lazy to write a matrix inversion routine ourselves. The installation of these libraries can be done by hand on desktop machines or using modules on HPC environments.
The recommended submit script uses Python and needs NetCDF4
#### Ubuntu/Debian Linux
Since Ubuntu 14.04LTS, the system libraries can be used:
```bash
sudo apt-get install gfortran libnetcdff5 libfftw3-dev libblas3 openmpi-bin gfortran python python-netcdf
```

#### MacOS X
We recommend using Homebrew to download and compile the necessary libraries. Be careful to add the *--with-fortran* argument while installing any library. Otherwise, only C libraries are installed. If you did install the libraries before without the *--with-fortran* argument, it may be necessary to remove and reinstall them.

#### HPC environments
Be sure to load modules for Fortran, NetCDF, FFTW and LAPACK.
On SuperMUC, the necessary commands are
```
module load fortran
module load netcdf/mpi
module load fftw
module load mkl
```

## AxiSEM wavefields
Using the code requires computation of a global seismic wavefield using AxiSEM. Please refer to the AxiSEM documentation on how to do this. An set of example wavefields with a dominant period of 40s can be downloaded:
```
wget https://www.geophysik.uni-muenchen.de/~staehler/kerner_wavefields.tar.bz2
tar -xvf kerner_wavefields.tar.bz2
```

## Installation
Change into the download directory and copy the included template files into the main directory
```bash
cd kerner
./copy_templates.sh
```
The file *make_mc_kerner.macros* allows you to modify the compiler name and compiler flags according to your system. To compile, use
```bash
make 
```
and afterwards 
``` bash
make check
```
to run a set of tests (wavefields downloaded in the previous step are necessary for the tests to complete).

## Usage
To run the code, a convenient submit script is available:
```bash
python submit.py -n NSLAVES JOB_NAME
```
This creates a run directory *JOB_NAME*, puts the default input files there and starts a MPI job with one master and NSLAVES worker tasks. In this simple form, it uses default values stored in the files *CMTSOLUTION*, *receivers.dat*, *filters.dat* and *stf_20s.dat* together with AxiSEM wavefields in ./wavefield/fwd and ./wavefield/bwd.  
Modifying the settings for the run is possible by setting one of the options of submit.py, which can be queried by
```bash
python submit.py -h
```

## Plot results
MC Kerner comes with a set of simple Python plotting scripts to visualize kernels and the seismograms they are based on. They can be found in the folder **pyfiles** and need a few extra programs to run:
 * Paraview (version >= 4.3)
 * Python packages matplotlib and wand

To run the scripts change into the job directory
```bash
cd JOB_DIR
```
and run
```bash
python ../pyfiles/plot_all_seismograms.py
pvpython ../pyfiles/plot_all_kernels.py
python ../pyfiles/create_composites.py
```
The last script creates a composite image like the one below:
![](https://www.geophysik.uni-muenchen.de/~staehler/kerner/composite_plot.png)




## FAQs
### Anaconda interference 
The python package manager Anaconda installs his own version of gfortran, which interferes with the system version and the MPI libraries. 
A solution is to replace the compiler name *FC* and the *MPIRUN* variable in *make_mc_kerner.macros* with absolute paths
```bash
FC = /usr/bin/mpif90
```
