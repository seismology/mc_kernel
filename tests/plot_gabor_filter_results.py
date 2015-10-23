
# coding: utf-8

# In[1]:

import filtering
import numpy as np
import matplotlib.pyplot as plt

# Load input time trace
signal_in_temp = np.loadtxt('../tests/gaborinput.txt')

len_signal = 1000
dt = 0.1
t = np.array(range(0, len_signal)) * dt

signal_in = np.zeros((len_signal, 1))
signal_in[:, 0] = signal_in_temp

# Define filter bank
pmax = (2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0)
tshift = (0.0, 2.5, 0.0, 2.5, 0.0, 2.5, 0.0, 5.0, 0.0, 5.0, 0.0, 5.0)
sigmaIfc = (0.3, 0.3, 0.5, 0.5, 0.7, 0.7, 0.3, 0.3, 0.5, 0.5, 0.7, 0.7)

fig = plt.figure(figsize=(15, 10))

# Load MC Kernel results for each filter and plot
for ifilt in range(0, len(pmax)):
    ax = fig.add_subplot(4, 3, ifilt + 1)
    signal_out_ref = filtering.gaborfilter(signal_in, dt=dt,
                                           pmax=pmax[ifilt],
                                           nscale=1, fmult=2,
                                           sigmaIfc=sigmaIfc[ifilt],
                                           npad=2048,
                                           tshift=(tshift[ifilt], 0))

    # Load MC Kernel results
    fnam = '../tests/gaboroutput_%02d.txt' % (ifilt+1)
    signal_out_kernel = np.loadtxt(fnam)
    diff = (signal_out_kernel[0:len_signal] -
            signal_out_ref[0][0:len_signal, 0, 0])

    # Plot Reference, MC Kernel and Difference * 1e6
    ax.plot(t, signal_out_ref[0][0:len_signal, 0, 0], 'r', label='PyFFProc')
    ax.plot(t, signal_out_kernel[0:len_signal], 'k', label='MC Kerner')
    ax.plot(t, diff*1e6, 'b', label='Difference x $10^6$')
    ax.set_xlim((45, 55))
    if (ifilt == 2):
        ax.legend()

fig.savefig('Gaborfilter_compare_test_reference.png', dpi=200)
