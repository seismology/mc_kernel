#!/usr/bin/env python
# coding: utf-8
"""
Test MC Kernel results against reference solution (for Travis CI)
Author: Simon Stähler, LMU München
"""
import os
import subprocess
from netCDF4 import Dataset
import numpy.testing as npt
import numpy as np
import glob
import argparse


helptext = 'Test MC Kernel results against reference solution'
parser = argparse.ArgumentParser(description=helptext)

helptext = "Run directory. The code will compare results in \n" + \
           "./tests/reference_solutions/RUNDIR with ./runs/RUNDIR"
parser.add_argument('rundir', help=helptext)

helptext = "Relative tolerance for kernel values"
parser.add_argument('--rtol', help=helptext, type=float, default=1e-5)

helptext = "Absolute tolerance for kernel values"
parser.add_argument('--atol', help=helptext, type=float, default=1e-15)

args = parser.parse_args()

ref_dirs = 'tests/reference_solutions'
dirname_ref = os.path.join(ref_dirs, args.rundir)
dirname_res = os.path.join('runs', args.rundir)

# Test ASCII files that should be absolutely identical
filelist = ['kerner_kernel.xdmf']

for file in filelist:
    filename = os.path.split(file)[-1]
    retcode = subprocess.call(
        ['diff',
         os.path.join(dirname_ref, filename),
         os.path.join(dirname_res, filename)])
    npt.assert_equal(retcode, 0,
                     err_msg='%s not equal to reference' % filename)

# Test seismograms, where numerical rounding errors might be present
filelist_smgr = glob.glob(os.path.join(dirname_ref, 'Seismograms', '*'))

for file in filelist_smgr:
    filename = os.path.split(file)[-1]
    seis_ref = np.loadtxt(os.path.join(dirname_ref, 'Seismograms', filename))
    seis_res = np.loadtxt(os.path.join(dirname_res, 'Seismograms', filename))
    npt.assert_allclose(seis_ref, seis_res, rtol=1e-8, atol=1e-16)

# Test kernel values, where rounding errors might be huge.
with Dataset(os.path.join(dirname_ref, 'kerner_kernel.nc'), 'r') as nc_ref:
    with Dataset(os.path.join(dirname_res, 'kerner_kernel.nc'), 'r') as nc_res:
        kernel_ref = np.array(nc_ref.groups['cell_data'].
                              variables['kernel'][:])
        kernel_res = np.array(nc_res.groups['cell_data'].
                              variables['kernel'][:])

        npt.assert_allclose(kernel_ref, kernel_res,
                            atol=args.atol, rtol=args.rtol)

        error_ref = nc_ref.groups['cell_data'].variables['error'][:]
        error_res = nc_res.groups['cell_data'].variables['error'][:]

        npt.assert_allclose(error_ref, error_res,
                            atol=args.atol, rtol=args.rtol)
