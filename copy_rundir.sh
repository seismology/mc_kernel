#!/bin/bash
echo "Creating directory $2 with same input files as in $1"
mkdir $2
cp $1/inparam_basic $1/CMTSOLUTION $1/filters.dat $1/receiver.dat $2



