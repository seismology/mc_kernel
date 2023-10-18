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

## Prerequisite
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
Using the code requires computation of a global seismic wavefield using AxiSEM. Please refer to the AxiSEM documentation on how to do this. A set of example wavefields with a dominant period of 40s can be downloaded:
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
This creates a run directory *JOB_NAME*, puts the default input files there and starts an MPI job with one master and NSLAVES worker tasks. In this simple form, it uses default values stored in the files *CMTSOLUTION*, *receivers.dat*, *filters.dat* and *stf_20s.dat* together with AxiSEM wavefields in ./wavefield/fwd and ./wavefield/bwd.  
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

Allowed values for the filter type are:
- `Gabor`: Log-frequency Gabor filter. Parameters: center period, width as ratio of center period, time shift in seconds (typically 0), *ignored*
- `Butterw_BP`: Bandpass Butterworth. Parameters: Lower cutoff period, Upper cutoff period, order of filter, *ignored*
- `Butterw_LP`: Lowpass Butterworth. Parameters: Upper cutoff period, *ignored*, order of filter, *ignored*
- `Butterw_HP`: Highpass Butterworth. Parameters: Lower cutoff period, *ignored*, order of filter, *ignored*
- `ident`: Do not filter at all

All Butterworth filters can be turned into a zero-phase filter by adding `_ZP` to the type. Note that this doubles the filter order.



## FAQs
### Anaconda interference 
The python package manager Anaconda installs his own version of gfortran, which interferes with the system version and the MPI libraries. 
A solution is to replace the compiler name *FC* and the *MPIRUN* variable in *make_mc_kerner.macros* with absolute paths
```bash
FC = /usr/bin/mpif90
```

## Create new AxiSEM wavefields
### Prerequisites
#### AxiSEM
Install the latest version of AxiSEM from github. The versions published on axisem.info are rather old.
### Some more comments
Remember: To compute a kernel, two groups with a total of ***6*** *AxiSEM* runs are necessary:
- 4 ***forward*** runs with double couple sources, for the 4 basis functions, into which an arbitrary wavefield of a moment tensor source in an axially symm
  axially symmetric model can be decomposed (MXX_P_MYY     MXY_MXX_M_MYY MXZ_MYZ       MZZ)
- 2 ***backward*** runs with single force sources, one vertical (PZ), one horizontal (PX).

AxiSEM contains a Python script to group the forward and the backward runs respectively. Note the following:

- The forward run needs to be at the ***depth of the quake***. This means that you need to do a new ***forward*** run 
if you want to investigate another quake, or general sensitivity of quakes in a different depth.
- The backward run needs to have its source at the depth of the receiver. This is 0 km in all plausible cases unless your
seismometer is in a really deep case (which it probably isn't). So typically, you will only do one ***backward*** run 
per model and frequency.
- After the AxiSEM run, the wavefield files need to be transformed, from contiguous in space (as created during the run) 
to contiguous in time (as needed by MC Kernel or Instaseis). This process is triggered by the submit script.
- The wavefields which AxiSEM stores can become really big (about 1 TB for a 1800 s long run at 2.5 second dominant 
period). The size scales linearly with duration of the simulation and cubic with maximum frequency. There is a variable 
in axisem/inparam_advanced which allows to restrict the wavefields to a certain depth range. This is typically **NOT** 
desired if you want to compute kernels!


### AxiSEM settings
go to your AxiSEM root directory ($AXISEM)

    cd $AXISEM

and copy the template input files

    ./copy_templates.sh

then modify two input files. 

#### Simulation length
In **inparam_basic** set the simulation length to about 100 seconds after the last kernel time
window that you want to use. If you want to look at 
 - S-waves (and SS, SKS) until antipodal distances: one hour or 3600s will be enough. This is the reasonable maximum.
 - P-waves and S-waves until 90 degree distance: about 40 minutes or 2400s. 
 - P-waves only until 90 degree distance: 1200s.
Since the AxiSEM runtime and wavefield storage size scales linearly with seismogram length, and MCKernel runtime scales
at least linearly with seismogram length, you can save a bit of time here.


    SIMULATION_LENGTH  1200

There are a number of other settings in inparam_basic, which at this stage should be left alone.
#### Wavefield output settings
In **inparam_advanced** check these parameters 

```
KERNEL_WAVEFIELDS  true
KERNEL_COLAT_MIN      0.
KERNEL_COLAT_MAX    180.
KERNEL_RMIN           0.
KERNEL_RMAX        7000.
```

KERNEL_RMIN should be set to 0 km (or the radius of the innermost cell in your kernel mesh), KERNEL_RMAX to at least the 
radius of your outermost mesh cell (typically the radius of the planet). It can be bigger, so 7000 is conservative for
the Earth radius in kilometers.
Again: Leave the other settings untouched.


### Run AxiSEM
Setting up the 6 runs is a bit involved, so there is a Python script which does most of the work for you:

````
python ./submit.py -h
usage: submit.py [-h] [--nrad NRAD] [--ntheta NTHETA] [-j JOB_TYPE] [-r RUN_TYPE] [--max_colat MAX_COLAT]
                 [--max_depth MAX_DEPTH] [--src_depth SRC_DEPTH] [--ncl NCL] [-w WALL_TIME] [-m MAIL_ADDRESS] [-a ACCOUNT]
                 job_name mesh_file mesh_period

Create AxiSEM run and submit job.

positional arguments:
  job_name              Job directory name. 
  mesh_file             Mesh file. Choose path to external mesh file or one
                        of the following AxiSEM internal models:
                        'prem_iso', 'prem_iso_solid', 'prem_iso_onecrust',
                        'prem_iso_light', 'prem_iso_solid_light',
                        'prem_ani', 'prem_ani_onecrust', 'prem_ani_light',
                        'ak135', 'ak135f', 'iasp91'
  mesh_period           Mesh period 

optional arguments:
  -h, --help            show this help message and exit
  --nrad NRAD           Number of radial slices (default: 1)
  --ntheta NTHETA       Number of theta slices (default: 2)
                        Set to 0 to use the maximum number possible for this
                        model and frequency
  -j JOB_TYPE, --job_type JOB_TYPE
                        Job type (local, Euler, CSCS Daint)
  -r RUN_TYPE, --run_type RUN_TYPE
                        Run type (Instaseis, MC Kernel)
                        Options: 
                           instaseis (default) 
                           mckernel_fwd
                           mckernel_bwd
  --src_depth SRC_DEPTH
                        Source depth in kilometer
  --ncl NCL             Number of coarsening layers (default: 1)
  -w WALL_TIME, --wall_time WALL_TIME
                        Wall time for the solver in hours
  -m MAIL_ADDRESS, --mail_address MAIL_ADDRESS
                        Mail adress for notifications
  -a ACCOUNT, --account ACCOUNT
                        Daint project accoun
````

To compute wavefields for an MC Kernel run, use the -r mckernel_fwd or -r mckernel_bwd run types. For example, to 
fire off runs for a Mesh with 20 second period in IASP91, that uses 8 cores (use less cores if your computer has less). 
```
python ./submit.py FWD_iasp91 iasp91 20 -r mckernel_fwd --ntheta 8
```
The first argument *FWD_iasp91* is a name for the run, the second argument chooses the model iasp91, the third sets 
the minimum period, the fourth makes it a forward run with 4 moment sources, the last chooses 8 cores.
To give you an idea: This run should take about ten minutes in total on a Macbook Pro 2019.
Monitor its progress with 

    tail -f ./runs/FWD_iasp91/MZZ/OUTPUT

After the MZZ run is done, check the output of MXX_P_MYY, MXY_MXX_M_MYY, MXZ_MYZ. After the 4 simulations are done, 
the wavefield database is automatically transformed to time-contiguous.

The backward run is triggered with
```
python ./submit.py BWD_iasp91 iasp91 20 -r mckernel_bwd --ntheta 8
```
and should take half as long.

### Copy the database files
After the backward run is done, check that there are database files created
```
ls -lh runs/FWD_iasp91/FWD_iasp91_database/MZZ/Data
```
This should show a file ordered_output.nc4 with roughly 100 MB in size.

Next, copy the databases to the MC Kernel wavefield directory
```
cp -r runs/FWD_iasp91/FWD_iasp91_database ../mc_kernel/wavefields
cp -r runs/BWD_iasp91/BWD_iasp91_database ../mc_kernel/wavefields
```

You can delete the AxiSEM run dir now, if you need to save diskspace. Otherwise, you're done and can run MC Kernel!
