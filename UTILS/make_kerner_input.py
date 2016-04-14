
# coding: utf-8

# In[1]:

from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
from obspy import read_events
from obspy import read_inventory
# import os
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs

helptext = 'Create MC Kernel receiver and source files.'
formatter_class = argparse.RawTextHelpFormatter
parser = argparse.ArgumentParser(description=helptext,
                                 formatter_class=formatter_class)

parser.add_argument('-e', '--event', help='Path to QuakeML file', required=True)
parser.add_argument('-r', '--receivers', nargs='+', required=True,
                    help='Path to StationXML file')
parser.add_argument('--nfilter', required=True,
                    help='Number of filters', type=int)
parser.add_argument('--f0', required=True,
                    help='Central frequency of shortest filter', type=float)

parser.add_argument('--phases', nargs='+', default='P',
                    help='Phases to calculate kernels for')

parser.add_argument('--parameters', nargs='+', default='P',
                    help='Model parameters to calculate kernels for')

args = parser.parse_args()

model_params = args.parameters
nparam = len(model_params)

# Read event file
events = read_events(args.event)
event = events[0].origins[0]
events.write('CMTSOLUTION', format='CMTSOLUTION')

# Read stations file
inventory = read_inventory(args.receivers[0])
for recfile in args.receivers[1:-1]:
    inventory += read_inventory(recfile)

filter_name = []
filter_length = []
filter_freq = []

for ifilter in range(0, args.nfilter):
    filter_name.append('band%02d' % ifilter)
    filter_freq.append(args.f0 * 2.0**(ifilter*0.5))
    filter_length.append(args.f0 * 2.0**(ifilter*0.5) * 1.5)

with open('filters.dat', 'w') as fid:
    fid.write('%d   Number of filters\n' % args.nfilter)
    for ifilter in range(0, args.nfilter):
        fid.write('%s %s %4.1f 0.5 0 0\n' % (filter_name[ifilter],
                                             'Gabor',
                                             filter_freq[ifilter]))

# Version to loop over stations in inventory
model = TauPyModel(model='ak135')

# Length of wavefield database
t_max = 1800.

stats = []
stat_names = []
# Create list of stations, independent of network
for network in inventory:
    for stat in network:
        # Remove duplicate stations
        if (stat.code not in stat_names):
            stats.append(stat)
            stat_names.append(stat.code)

# Sort stations by distance to source
stats_sort = sorted(stats, key=lambda stat: locations2degrees(stat.latitude,
                                                              stat.longitude,
                                                              event.latitude,
                                                              event.longitude))
nrec = len(stats_sort)

with open('receiver.dat', 'w') as fid:
    fid.write('%d   # nr of receivers\n' % nrec)
    fid.write('Z    # component\n')

    phase_list = args.phases #['P', 'Pdiff', 'PP', 'S', 'Sdiff', 'SS']
    # 'PPP', 'SSS'] #, 'PKiKP' ] #, 'PcP', 'PPP', 'PKiKP', 'PKIKP', 'PKP']

    for stat in stats_sort:
        distance = locations2degrees(stat.latitude, stat.longitude,
                                     event.latitude, event.longitude)
        # print distance
        tt = model.get_travel_times(distance_in_degree=distance,
                                    source_depth_in_km=event.depth*1e-3,
                                    phase_list=phase_list)

        # Count number of phases for this distance
        nphase = 0
        phases = []
        iphases = []
        for i in range(0, len(tt)):
            if tt[i].name in phase_list:
                if (not tt[i].name in phases) and \
                   (tt[i].time + max(filter_length) < t_max):
                    iphases.append(i)

                phases.append(tt[i].name)

        nphase = len(iphases)

        fid.write('%s  %8.3f %8.3f %d\n' % (stat.code, stat.latitude,
                                            stat.longitude,
                                            nphase * args.nfilter * nparam))
        print('%6s, %8.3f degree, %d Kernel\n' % (stat.code, distance,
                                                  nphase * args.nfilter * nparam))
        for i in iphases:
            if tt[i].name in phase_list:
                travel_time = tt[i].time
                for iband in range(0, args.nfilter):
                    for param in model_params:
                        fid.write('%s_%s_%02d %s CC %8.3f %8.3f %s\n' %
                                (tt[i].name,
                                 param,
                                 iband+1,
                                 filter_name[iband],
                                 travel_time - 5,
                                 travel_time + filter_length[iband]-5,
                                 param) )

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(111, projection=ccrs.Robinson())
# make the map global rather than have it zoom in to
# the extents of any plotted data
ax.set_global()

#ax.stock_img()
ax.coastlines()

ax.plot(event.longitude, event.latitude, '*', color='red',
        markersize=15, transform=ccrs.PlateCarree())

# add a circle
for dist in (30, 90, 120):
    circle = mpatches.Circle((event.longitude, event.latitude), dist,
                              ec="black", fc="None", ls="dotted", transform=ccrs.PlateCarree())
    ax.add_patch(circle)


stats = []
stat_names = []
# Create list of stations, independent of network
for network in inventory:
    #print(network)
    for stat in network:
        if (not stat.code in stat_names):
            stats.append(stat)
            stat_names.append(stat.code)

# Sort stations by distance to source
stats_sort = sorted(stats[0:nrec], key=lambda stat:locations2degrees(stat.latitude, stat.longitude,
                                                                     event.latitude, event.longitude))
nrec = len(stats_sort)

for stat in stats_sort:
    ax.plot(stat.longitude, stat.latitude, 'v', color='darkgreen',
            transform=ccrs.PlateCarree(), markersize=10)
ax.legend(fontsize=12, numpoints=1, markerscale=2.0)

fname = './stations.png'
fig.savefig(filename=fname, 
            bbox_inches='tight',
            dpi=300)
#
#




# import numpy as np
# lat1 = stat.latitude * np.pi / 180
# lon1 = stat.longitude * np.pi / 180
# lat2 = event.latitude * np.pi / 180
# lon2 = event.longitude * np.pi / 180
# v1 = (np.cos(lon1)*np.cos(lat1), np.sin(lon1)*np.cos(lat1), np.cos(lat1))
# v2 = (np.cos(lon2)*np.cos(lat2), np.sin(lon2)*np.cos(lat2), np.cos(lat2))
#
# print np.cross(v1, v2)

