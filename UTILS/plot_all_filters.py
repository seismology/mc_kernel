#!/usr/bin/env python
# coding: utf-8
"""
Plot all filters in this rundir
Author: Simon Stähler, LMU München
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import glob

filt_dir = './Filters/'
filt_plot_dir = './Filter_plots'

filt_list = glob.glob(os.path.join(filt_dir, 'filterresponse_stf_*'))

os.mkdir(filt_plot_dir)

for filename in filt_list:
    filt_split = os.path.split(filename)[-1].split('_')
    filt_name = filt_split[2]
    freq_1 = float(filt_split[3])
    freq_2 = float(filt_split[4])

    print('Filter: %s (%f, %f)' % (filt_name, freq_1, freq_2))

    filt_data = np.loadtxt(filename)

    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    ax.plot(1. / filt_data[:, 0],
            np.sqrt(filt_data[:, 1]**2 + filt_data[:, 2]**2),
            'k', lw=3, zorder=10, label='Filterresponse')
    ax.plot(1. / filt_data[:, 0],
            np.sqrt(filt_data[:, 3]**2 + filt_data[:, 4]**2),
            'k', linestyle='dashed', lw=3, zorder=11,
            label='Filterresponse, including STF')

    # Plot all other filters in grey thin lines
    for filename in filt_list:
        filt_data = np.loadtxt(filename)
        ax.plot(1. / filt_data[:, 0],
                np.sqrt(filt_data[:, 1]**2 + filt_data[:, 2]**2),
                color='lightgrey', lw=1, zorder=0)
        ax.plot(1. / filt_data[:, 0],
                np.sqrt(filt_data[:, 3]**2 + filt_data[:, 4]**2),
                color='lightgrey', linestyle='dashed', lw=1, zorder=1)

    ax.set_xlim((0, 500))
    ax.set_ylim((0, 1))

    ax.set_xlabel('Period / s')
    fig.savefig(os.path.join(filt_plot_dir,
                             '%s_%6.4f_%6.4f.png' % (filt_name,
                                                     freq_1, freq_2)))

    plt.close(fig)
