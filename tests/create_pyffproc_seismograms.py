
# coding: utf-8

#
# simple and convoluted script creating instaseis seismograms,
# filtering them with the pyffproc loggabor routine and storing
# the filtered seismogram for comparison with mckernel
#

import instaseis
import numpy as np
from filtering import bandpassfilter
import matplotlib.pyplot as plt

plot = True

# Create STF:
dt = 1.0
t_half = 10.0
t = np.arange(0,100,dt)
stf = 2 * np.exp(-(t/t_half)**2) * t / t_half**2
nsample = len(t)

f = open('./stf_pyffproc.dat', 'w')

f.write('%d %f\n'%(nsample, dt))
for y in stf:
    f.write('%f\n'%(y))

f.close()

# Open instaseis db
db = instaseis.open_db("/scratch/auerl/wavefield_pyffproc_10s_647km_prem_ani/bwd/")
#db = instaseis.open_db('../wavefield/bwd')

receiver = instaseis.Receiver(latitude=-51.68, longitude=-58.06, network="AB", station="EFI")
source = instaseis.Source(latitude   = -13.8200,
                          longitude  = -67.2500,
                          depth_in_m = 647100,
                          m_rr = -7.590000e+27 / 1E7,
                          m_tt =  7.750000e+27 / 1E7,
                          m_pp = -1.600000e+26 / 1E7,
                          m_rt = -2.503000e+28 / 1E7,
                          m_rp =  4.200000e+26 / 1E7,
                          m_tp = -2.480000e+27 / 1E7,
                          time_shift=None, 
                          sliprate=stf, 
                          dt=dt)
source.resample_sliprate(db.info.dt, db.info.npts)

# create velocity and displacement seismograms
st = db.get_seismograms(source, receiver, reconvolve_stf=True,
                        kind="displacement", remove_source_shift=False)
st.filter('lowpass', freq=1./10.0)
tr_disp = st[0]
st = db.get_seismograms(source, receiver,  reconvolve_stf=True,
                        kind="velocity", remove_source_shift=False)
st.filter('lowpass', freq=1./10.0)
tr_velo = st[0]

# bandpass the seismograms
len_orig=np.size(tr_disp.data)
tr_velo_bandpassed = bandpassfilter(tr_velo, 'log-gabor', 1024, 8, 30, 1.4142, 0.5)
tr_disp_bandpassed = bandpassfilter(tr_disp, 'log-gabor', 1024, 8, 30, 1.4142, 0.5)

if plot:
    f, axarr = plt.subplots(3, sharex=True,sharey=True)

# plot and export
for i in np.arange(0,4,1):

    f = open('seism_ref_raw_EFI_P_0'+str(i+1), 'w')
    f.write('%d\n'%len_orig)
    for y in tr_disp.data:
        f.write('%e\n'%y)
    f.close()
    
    f = open('seism_ref_EFI_P_0'+str(i+1), 'w')
    f.write('%d\n'%len_orig)
    for y in np.arange(0,len_orig):
        a=tr_disp_bandpassed[0][y, 0, i]
        b=tr_velo_bandpassed[0][y, 0, i]        
        f.write('%e %e\n' % (a,b))
    f.close()

    if plot:
        c=np.genfromtxt('./Seismograms/seism_EFI_P_0'+str(i+1))
        time=c[:,0]
        axarr[0].plot(time,tr_velo_bandpassed[0][0:len_orig, 0, i])
        axarr[1].plot(time,c[:,2])
        axarr[2].plot(time,tr_velo.data)
        axarr[0].set_title('Gabor filter results in pyffproc')
        axarr[1].set_title('Gabor filter results in mckernel')
        axarr[2].set_title('Raw unfiltered seismogram')

if plot:
    plt.show()    
