
# coding: utf-8



import instaseis
import numpy as np


# Create STF:
dt = 1.0
t_half = 20.0
t = np.arange(0,100,dt)
stf = 2 * np.exp(-(t/t_half)**2) * t / t_half**2
nsample = len(t)

f = open('../stf_test.txt', 'w')

f.write('%d %f\n'%(nsample, dt))
for y in stf:
    f.write('%f\n'%(y))

f.close()

print sum(stf*dt)



db = instaseis.open_db('../wavefield/bwd')
receiver = instaseis.Receiver(latitude = 30, longitude = 0, network='MC', station='kerner')
source = instaseis.Source(latitude = 90.0, 
                          longitude = 0.0, 
                          depth_in_m=0.0, 
                          m_rr=1.0e13,
                          m_tt=1.0e13, 
                          m_pp=1.0e13, 
                          m_rt=0.0, 
                          m_rp=0.0, 
                          m_tp=0.0,
                          time_shift=None, 
                          sliprate=stf, 
                          dt=dt)
source.resample_sliprate(db.info.dt, db.info.npts)



st = db.get_seismograms(source, receiver, reconvolve_stf=True, kind="displacement")
st.filter('lowpass', freq=1./40.0)
st.plot()
tr_Z = st[0]

tr_Z.data *= 1e20

f = open('seismogram_ref_T001_P_BW', 'w')

f.write('%d\n'%tr_Z.stats.npts)

for y in tr_Z.data:
    f.write('%e\n'%y)

f.close()

