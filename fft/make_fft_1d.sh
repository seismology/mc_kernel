#!/bin/bash

#gfortran -Wall -fbounds-check -pedantic -fbacktrace -c -I $HOME/local/include/ -I /usr/include fft_1D.f90
gfortran -Wall -fbounds-check -pedantic -fbacktrace -c -I /usr/include fft_1D.f90
gfortran fft_1D.o -o fft_1D  -lfftw3 -lfftw3f -Wall -pedantic -fbacktrace
