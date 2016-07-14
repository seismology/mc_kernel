#!/usr/bin/env python
# coding: utf-8

## Update kerner time windows
# Reads all time windows in receiver.dat, plots them with instaseis and allows
# the user to update them by clicking into a matplotlib window.
# Could use some more functionality, like plotting vlines after the user clicks.

### Initialize instaseis and matplotlib

# In[1]:

import instaseis
import matplotlib.pyplot as plt
import numpy as np
import os

kerner_dir = '/home/staehler/kerner/'


### Find directory of backward run (for instaseis)

# In[2]:


f = open(os.path.join(kerner_dir, 'inparam_basic'), 'r')
str = f.readline()
#print str
dir_not_found_yet = True
while dir_not_found_yet and str!='':
    try:
        if str.index('BWD_DIR') == 0:
            dirname = str.split()[1]
            #Remove hyphens, if there
            dirname = dirname.lstrip("'")
            dirname = dirname.rstrip("'")
            print dirname
            f.close()
            dir_not_found_yet = False
    except ValueError:
        str = f.readline()




### Load instaseis database

# In[3]:

db = instaseis.open_db(os.path.join(kerner_dir, dirname))
print db


### Load CMTSOLUTION file

# In[4]:

src = instaseis.Source.parse(os.path.join(kerner_dir, 'CMTSOLUTION'))
print src


### Parse receivers.dat file

# In[ ]:

def onclick(event):
    if event.button==3:
        print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            event.button, event.x, event.y, event.xdata, event.ydata)
        new_time_window.append(event.xdata)

        ax.vlines(event.xdata, -1, 1)
        #if len(new_time_window) == 2:
        #    fig.canvas.mpl_disconnect(cid)
    return new_time_window

#
f = open(os.path.join(kerner_dir, 'receiver.dat'))
f_new = open(os.path.join(kerner_dir, 'receiver_new.dat'), 'w')

# Read number of receivers
str_line = f.readline()
f_new.write(str_line)
nrec = int(str_line.split()[0])
print 'Number of receivers: %d'%nrec

# Read seismogram component
str_line = f.readline()
seis_cmp = str_line.split()[0]
f_new.write(str_line)
print 'Seismogram component: %s'%seis_cmp

for irec in range(0, nrec):
    str_line = f.readline()
    rec_name = str_line.split()[0]
    rec_lat  = float(str_line.split()[1])
    rec_lon  = float(str_line.split()[2])
    nkernel  = int(str_line.split()[3])
    f_new.write(str_line)

    print 'Receiver: %s, coordinates: (%f, %f), %d kernels'%(rec_name, rec_lat, rec_lon, nkernel)
    rec = instaseis.Receiver(latitude=rec_lat, longitude=rec_lon)
    st = db.get_seismograms(source=src, receiver=rec, components=seis_cmp)
    tr = st[0]

    for ikernel in range(0, nkernel):
        str_line = f.readline()
        kernel_name = str_line.split()[0]
        filter_name = str_line.split()[1]
        misfit_name = str_line.split()[2]
        time_window_start = float(str_line.split()[3])
        time_window_stop  = float(str_line.split()[4])
        model_param = str_line.split()[5]

        print 'Kernel %s, filter %s, time window (%8.2fs %8.2fs)'%(kernel_name, filter_name, time_window_start, time_window_stop)
        tr_time_window = st[0].slice(starttime=tr.stats.starttime + time_window_start, endtime=tr.stats.starttime + time_window_stop)
        time_window_times = tr_time_window.times() + float(str_line.split()[3])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(tr.times(), tr.data)
        ax.plot(time_window_times, tr_time_window.data, 'r')

        ax.set_xlabel('time / seconds')
        ax.set_ylabel('amplitude')
        ax.set_title('Kernel %s, time window (%8.2fs %8.2fs)'%(kernel_name, time_window_start, time_window_stop))
        cid = fig.canvas.mpl_connect('button_press_event', onclick)

        print 'Mark new time window with two right mouse button clicks. Close window, when done'
        new_time_window = []

        plt.show()
        new_time_window = np.array(new_time_window)



        if len(new_time_window)==2:
            str_output = 'Do you want to update the time windows for this kernel to (%8.2fs %8.2fs) [y/n]?'%(new_time_window[0], new_time_window[1])
            update_time_window = raw_input(str_output)
            if update_time_window=='y':
                time_window_start = new_time_window[0]
                time_window_stop  = new_time_window[1]

                new_line = '%s %s %s %10.4f %10.4f %s\n'%(kernel_name, filter_name, misfit_name,
                                                          time_window_start, time_window_stop, model_param)

        f_new.write(new_line)


f.close()
f_new.close()



