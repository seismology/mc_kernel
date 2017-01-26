#!/bin/sh
# runtests.sh --
#    Bourne shell script to control a program that uses funit
#    Name of the program: first argument
#
#    $Id: runtests.sh,v 1.2 2008/01/26 11:15:10 arjenmarkus Exp $
#
set -e 
#echo "Compiling code"
#make -s
#
#
##if [ ! -d test_wavefields ] ; then
#    echo "Downloading test wavefield files"
#    wget http://geophysik.uni-muenchen.de/~staehler/test_wavefields.tar.gz
#    echo "Unpacking test wavefield files"
#    tar -xf test_wavefields.tar.gz
#fi
#
rm -f netcdf_out_*.tmp; 
rm -f OUTPUT_test
rm -f mckernel_tests.log
rm -f fail_messages.txt
rm -rf Seismograms
mkdir Seismograms
rm -rf Filters
mkdir Filters


echo "Running test"
set +e 
echo ALL >ftnunit.run
rm -f errors

chk=1
until test ! -f ftnunit.lst -a $chk -eq 0 ; do
    chk=0
    #valgrind --tool=memcheck $1 $2 $3 $4 $5 $6 $7 $8 $9 >>runtests.log 2>&1
    mpirun -n 1 $1 $2 >>mckernel_tests.log 2>&1
done
set -e 

# Present some test results
echo "##################################################"
tail -n 5 mckernel_tests.log
if grep -C 1 "assertion failed" ./mckernel_tests.log > fail_messages.txt
then
  echo 
  echo "##################################################"
  echo "FAILURES: At least one test failed"
  echo "Details of failed tests:"
  cat fail_messages.txt
  echo 
  echo "See tests/mckernel_tests.log for details"
  echo "##################################################"
  exit 1
elif [ -f errors ]
then
  echo 
  echo "##################################################"
  echo "ERRORS: At least one test crashed"
  echo "See tests/mckernel_tests.log for details"
  echo "##################################################"
  exit 1
else
  echo 
  echo "##################################################"
  echo "All tests successful"
  echo "##################################################"
fi
