#!/bin/bash

gfortran -Wall -pedantic -fbacktrace -c -I $HOME/local/include/ -I /usr/include fft_type.f90
gfortran fft_type.o -o fft_type -L $HOME/local/lib -Wl,-rpath,$HOME/local/lib -L /usr/lib -lm -lfftw3 -lfftw3f -Wall -pedantic -fbacktrace
