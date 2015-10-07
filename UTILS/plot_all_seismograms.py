# coding: utf-8
"""
Plot all seismograms in this rundir with a zoom window on the
respective time windows
Author: Simon Stähler, LMU München
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import glob

seis_dir = './Seismograms/'
seis_plot_dir = './Seismogram_plots'

smgr_list = glob.glob(os.path.join(seis_dir, 'seism_*'))

os.mkdir(seis_plot_dir)

for filename in smgr_list:
    kernel = os.path.split(filename)[1][6:]
    if kernel[0:3] in ('cut', 'raw'):
        continue

    print kernel

    smgr_name = os.path.join(seis_dir, 'seism_%s' % kernel)
    cut_name = os.path.join(seis_dir, 'seism_cut_%s' % kernel)

    smgr = np.loadtxt(smgr_name)
    cut = np.loadtxt(cut_name)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    ax.plot(smgr[:, 0], smgr[:, 1], 'k', lw=3)
    ax.plot(cut[:, 0], cut[:, 1], 'r', lw=3)
    norm = max(abs(smgr[:, 1]))
    ax.set_xlim((0, 2400))
    ax.set_ylim(-norm*1.5, norm*1.0)

    ax.set_xlabel('Time')

    # this is an inset axes over the main axes
    axins = fig.add_axes([.2, .2, .3, .3], axisbg='y')
    axins.plot(smgr[:, 0], smgr[:, 1], 'k', lw=3)
    axins.plot(cut[:, 0], cut[:, 1], 'r', lw=3)
    axins.set_xlim((cut[0, 0]-10, cut[-1, 0]+10))
    norm = max(abs(cut[:, 1]))
    axins.set_ylim(-norm*1.2, norm*1.2)

    axins.set_title('Zoom to phase')

    fig.savefig(os.path.join(seis_plot_dir, '%s.png' % kernel))
    plt.close(fig)
