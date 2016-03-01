# coding: utf-8

import numpy as np
import os
import subprocess
import glob
import matplotlib.pyplot as plt


def tail(f, n, offset=0):
    cmd = ["tail", "-n", str(n), f]
    tail_output = subprocess.check_output(cmd)
    tail_output.split('\n')

    return lines

dirnam = '.'
output_files = glob.glob(os.path.join(dirnam, './OUTPUT*'))
output_files.sort()

fnam = output_files[1]

cmd = ["grep", "CLOCKS", fnam]
tail_output = subprocess.check_output(cmd, universal_newlines=True)
lines = str(tail_output).split("\n")[1:-2]
ntimers = len(lines)
timer_name = []

# Get the names of the timers
for line in lines:
    start = 1
    # Remove the '-' that some timer names have
    for i in range(1, 3):
        if line.split()[i] == '-':
            start += 1
    timer_name.append(" ".join(line.split()[start:-4]))

t_per_call = np.zeros(ntimers)
t_total = np.zeros(ntimers)
ncalls = np.zeros(ntimers)

# Go through each slave Output file and get the timing information
for fnam in output_files[1:-1]:
    cmd = ["grep", "CLOCKS", fnam]
    tail_output = subprocess.check_output(cmd, universal_newlines=True)

    lines = str(tail_output).split("\n")[1:-2]
    if len(lines)>0:
      for iline in range(0, ntimers):
          line = lines[iline]
          timing = line.split()[-4:-1]
          ncalls[iline] += int(timing[0])
          t_per_call[iline] += float(timing[1])
          t_total[iline] += float(timing[2])

interesting_timers = np.array([2, 9, 11, 12, 13, 15, 17])

timer_interesting = []
t_total_interesting = []
t_rest = 0.0
for i in range(0, ntimers):
    if i in interesting_timers:
        timer_interesting.append('%s, \n %8.1f CPUh' % (timer_name[i],
                                                        t_total[i] / 3600.))
        t_total_interesting.append(t_total[i])
    else:
        t_rest += t_total[i]

#t_total_interesting.append(t_rest)
#timer_interesting.append('%s, \n %8.1f CPUh' % ('Other',
#                                                t_rest / 3600.))

# Create pie chart
# IO colors are reddish, Computation blueish, MPI yellowish
colors = ['dodgerblue',    # FFT
          'firebrick',     # NetCDF
          'tomato',        # Buffer
          'lightskyblue',  # Calc_strain
          'aqua',          # Lagrange
          'cadetblue',     # Filtering
          'darkcyan',      # Integration
          'yellowgreen',   # MPI
          'grey']          # Rest
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111)
ax.set_aspect(1)
ax.pie(t_total_interesting,
       labels=timer_interesting,
       colors=colors, autopct='%1.1f%%', shadow=True)
ax.set_title('Total calculation cost:%8.1f CPUh' % \
             (np.sum(np.array(t_total_interesting)) / 3600.))
fig.savefig('Timing.pdf')
