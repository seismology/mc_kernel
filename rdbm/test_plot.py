#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np

#nstat = 1
#rdbm_dat = np.loadtxt('mzz_reci_test_fwd')
#
#for i in np.arange(nstat):
#    plt.plot(rdbm_dat[:,0], rdbm_dat[:,i+1], color='k', alpha=1)
#
#
#nstat = 1
#rdbm_dat = np.loadtxt('mzz_reci_test_bwd')
#
#for i in np.arange(nstat):
#    plt.plot(rdbm_dat[:,0], rdbm_dat[:,i+1] * .55, color='r', alpha=1)


nstat = 25

for i in np.arange(nstat):
    ref_dat = \
        np.loadtxt('/home/ex/local/src/axisem/SOLVER/50s_bwd_5gll_straintrace/Data_Postprocessing/SEISMOGRAMS/recfile_%04d_disp_post_mij_conv0000_Z.dat'
                    % (i+1,))
    rdbm_dat = np.loadtxt('seis_%03d' % (i+1,))

    maxi = np.max(ref_dat[:,1])
    plt.plot(ref_dat[:,0], ref_dat[:,1] / maxi + i, color='r', alpha=1)
    plt.plot(rdbm_dat[:,0]-147.5, rdbm_dat[:,1] / maxi + i, color='k', alpha=1)

plt.show()
