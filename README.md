# MC Kernel
[![DOI](https://zenodo.org/badge/21468/seismology/mc_kernel.svg)](https://zenodo.org/badge/latestdoi/21468/seismology/mc_kernel)
[![Build Status](https://travis-ci.org/seismology/mc_kernel.svg?branch=master)](https://travis-ci.org/seismology/mc_kernel)
[![Coverage Status](https://coveralls.io/repos/github/seismology/mc_kernel/badge.svg?branch=master)](https://coveralls.io/github/seismology/mc_kernel?branch=master)

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
MC Kernel needs a recent Fortran compiler (gfortran 4.8 or later), an MPI installation, the NetCDF library for accessing the wavefield files, the FFTW library for Fourier transforms. The installation of these libraries can be done by hand on desktop machines or using modules on HPC environments.
The recommended submit script uses Python and needs NetCDF4
#### Ubuntu/Debian Linux
Since Ubuntu 14.04LTS, the system libraries can be used:
```bash
- sudo apt-get install gfortran libnetcdff5 libnetcdf-dev libfftw3-dev openmpi-bin libopenmpi-dev
```

#### MacOS X
We recommend using Homebrew to download and compile the necessary libraries. Be careful to add the *--with-fortran* argument while installing any library. Otherwise, only C libraries are installed. If you did install the libraries before without the *--with-fortran* argument, it may be necessary to remove and reinstall them.

#### HPC environments
Be sure to load modules for Fortran, NetCDF, FFTW.
On SuperMUC, the necessary commands are
```
module load fortran
module load netcdf/mpi
module load fftw
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

## Input files
There are 3 types of inputs needed to compute kernels:
- Receivers (locations of seismograms)
- Filters (describing the frequency content)
- Kernels (given a receiver and a filter, what other parameters to use for the kernel)

Each Receiver is described by 3 parameters:
- **Name**: Descriptive name of the receiver (e.g. SEED snippet like CH.DAVOX or really anything)
- **Latitude**
- **Longitude**
- **Source time function**
  
Each kernel is described by 6 parameters
- **Receiver**: Each kernel is estimated for a seismogram recorded at one given receiver.
- **Kernel name**: Just a descriptive name
- **Misfit type**: Can be cross-correlation travel time (*CC*) or amplitude (*AM*)
- **Time window**: Beginning and End of time window in which the misfit is computed
- **Filter**: This contains the frequency content of the measurement.
- **Model parameter**: Physical parameter to which the kernel belongs (*vp*, *vs*, etc)

These inputs are contained in two files: *receiver.dat* and *filter.dat*
### receiver.dat
Example file describing a total of 4 kernels at 2 receivers; one in 30 and one in 60 degree distance. 
```
2                                      # Nr of receivers
Z                                      # seismogram component
30deg  60.0000   00.0000  src1 2       # Name, lat, lon, stf name, nr of kernels
P_40s Gabor_40 CC 345.0 405.0  vp      # Kernel name, Filter name, misfit, time window begin and end, model parameter
P_30s Gabor_30 CC 345.0 390.0  vp
60deg  30.0000   00.0000  src1 2       
P_40s Gabor_40 CC 580.0 640.0  vp      
P_30s Gabor_30 CC 580.0 625.0  vp
```
Detailed breakdown of each line:
```
2                                      # Nr of receivers
```
How many receiver blocks will be read. If there are more receiver blocks below than this number, the remaining blocks are all ignored.

```
Z                                      # seismogram component
```
Which seismogram component do we compute the kernel on. Note that M.C. Kernel needs to be run again for a different seismogram component.

```
30deg  60.0000   00.0000  src1 2       # Name, lat, lon, stf name, nr of kernels
P_40s Gabor_40 CC 345.0 405.0  vp      # Kernel name, Filter name, misfit, time window begin and end, model parameter
P_30s Gabor_30 CC 345.0 390.0  vp
```
Receiver block. Starts with one line defining the receiver location and the number of kernels described on it. Followed by one line per kernel.
Each kernel line contains the 6 entries described in the comment above. 
Allowed values for the misfit are:
- `CC`: Cross-correlation travel time misfit as defined in Sigloch & Nolet (2006)
- `AM`: Amplitude misfit as defined in Sigloch & Nolet (2006)

Allowed values for the model parameter are:
- `vp`: P-wave speed
- `vs`: S-wave speed
- `rho`: Density
- `vph`: P-wave speed in horizontal direction
- `vpv`: P-wave speed in vertical direction
- `vsh`: S-wave speed in horizontal direction
- `vsv`: S-wave speed in vertical direction
- `eta`: anisotropic parameter eta
- `phi`: anisotropic parameter phi
- `xi `: anisotropic parameter xi
- `lam`: First Lamé parameter
- `mu `: Second Lamé parameter

### filters.dat
Example filters.dat file describing 3 filters:
- Gabor filter at 40 sec period
- low-pass 2nd order Butterworth (at 20sec period)
- bandpass 6th order Butterworth (between 20 and 40 sec)

```
3
filter1 Gabor 40.0 0.5 0.0 0.0  
filter2 Butterw_HP 20.   0.0 2.0 0.0  
filter3 Butterw_BP 20.0 40.0 6.0 0.0  
```

Allowed values for the filter type are: `Butterw_BP` `Butterw_HP`,  `Butterw_LP`



## FAQs
### Anaconda interference 
The python package manager Anaconda installs his own version of gfortran, which interferes with the system version and the MPI libraries. 
A solution is to replace the compiler name *FC* and the *MPIRUN* variable in *make_mc_kerner.macros* with absolute paths
```bash
FC = /usr/bin/mpif90
```
