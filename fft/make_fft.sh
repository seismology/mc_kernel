#!/bin/bash

gfortran -Wall -pedantic -fbacktrace -c -I $HOME/local/include/ -I /usr/include fft.f90
gfortran fft.o -o fft -L $HOME/local/lib -Wl,-rpath,$HOME/local/lib -L /usr/lib -lm -lfftw3 -lfftw3f
