#!/bin/bash
#
# adapted for openmpi from
# https://www.lrz.de/services/compute/supermuc/tuning/gprof/

PDIR=profdir.$OMPI_COMM_WORLD_RANK
mkdir -p $PDIR
export GMON_OUT_PREFIX=$PDIR/gmon.out

./kerner >& "$@" 
