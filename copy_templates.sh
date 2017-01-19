#!/bin/bash
set -e

[ -f make_mc_kernel.macros ] && cp make_mc_kernel.macros make_mc_kernel.macros.OLD
if [ $1 == 'release' ]; then
  echo "Copying files for fast code"
  cp TEMPLATES/make_mc_kernel.macros_release.TEMPLATE make_mc_kernel.macros
else
  echo "Copying files for debugging"
  cp TEMPLATES/make_mc_kernel.macros.TEMPLATE make_mc_kernel.macros
fi

[ -f Makefile ] && cp Makefile Makefile.OLD
cp TEMPLATES/Makefile.TEMPLATE Makefile

[ -f CMTSOLUTION ] && cp CMTSOLUTION CMTSOLUTION.OLD
cp TEMPLATES/CMTSOLUTION.TEMPLATE CMTSOLUTION

[ -f inparam_basic ] && cp inparam_basic inparam_basic.OLD
cp TEMPLATES/inparam_basic.TEMPLATE inparam_basic

[ -f receiver.dat ] && cp receiver.dat receiver.dat.OLD
cp TEMPLATES/receiver.dat.TEMPLATE receiver.dat

[ -f filters.dat ] && cp filters.dat filters.dat.OLD
cp TEMPLATES/filters.dat.TEMPLATE filters.dat

[ -f stf_20s.dat ] && cp stf_20s.dat stf_20s.dat.OLD
cp TEMPLATES/stf_20s.dat.TEMPLATE stf_20s.dat

[ -f src/Makefile ] && cp src/Makefile src/Makefile.OLD
cp TEMPLATES/Makefile_src.TEMPLATE src/Makefile

