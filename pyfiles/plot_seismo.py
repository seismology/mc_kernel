"""
Simple script to plot two waveforms over each other.
Normal usage: 
look at the whole smgr and the selected part.
Example:
python plot_seismo.py seismogram_60deg_P_60 seismogram_cut_60deg_P_60
"""

import matplotlib.pyplot as plt
import numpy as np
import sys

if not sys.argv[1]:
    sys.exit("Usage: python plot_seismo.py add1 add2")
else:
    smgr = np.loadtxt(sys.argv[1])

if not sys.argv[2]:
    sys.exit("Usage: python plot_seismo.py add1 add2")
else:
    cut = np.loadtxt(sys.argv[2])

plt.ion()

plt.plot(smgr[:,0], smgr[:,1], 'b', lw=3)
plt.plot(cut[:,0], cut[:,1], 'r', lw=3)

plt.xlabel('Time', size=24, weight='bold')
plt.xticks(size=18, weight='bold')
plt.yticks(size=18, weight='bold')

plt.show()

raw_input('Press enter to exit...')
