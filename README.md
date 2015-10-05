![](https://raw.githubusercontent.com/sstaehler/kerner/master/doc/logo.png?token=AE0eeLqaWCLbz3zNLBkZjTzREWjCq0l0ks5WFu-hwA%3D%3D)

##Authors:
Simon Stähler, Martin van Driel, Ludwig Auer, Kasra Hosseini, Tarje Nissen-Meyer

## Download
To download the current developer version, use git clone
```
git clone https://github.com/sstaehler/kerner
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
make test 
```
to run a set of tests.

## AxiSEM wavefields
Using the code requires computation of a global seismic wavefield using AxiSEM. Please refer to the AxiSEM documentation on how to do this. An set of example wavefields with a dominant period of 40s can be downloaded:
```
wget https://www.geophysik.uni-muenchen.de/~staehler/kernel_wavefields_40s.tar.gz
tar -xvf kernel_wavefields_40s.tar.gz
```

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
![](https://raw.githubusercontent.com/sstaehler/kerner/master/doc/composite_plot.png?token=AE0eeJLSdHQBHr3dA-2R7myVym7q9jcYks5WFu-8wA%3D%3D)




## FAQs
### Anaconda interference 
The python package manager Anaconda installs his own version of gfortran, which interferes with the system version and the MPI libraries. 
A solution is to replace the compiler name *FC* and the *MPIRUN* variable in *make_mc_kerner.macros* with absolute paths
```bash
FC = /usr/bin/mpif90
```
