#!/usr/bin/env python
# coding: utf-8

import obspy
from obspy.taup import TauPyModel
from obspy.geodetics import locations2degrees
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cartopy.crs as ccrs

helptext = 'Create MC Kernel receiver and source files.'
formatter_class = argparse.RawTextHelpFormatter
parser = argparse.ArgumentParser(description=helptext,
                                 formatter_class=formatter_class)

parser.add_argument('-e', '--event', help='Path to QuakeML file')
parser.add_argument('-r', '--receivers', nargs='+',
                    help='Path to StationXML file')
parser.add_argument('--nfilter', required=True,
                    help='Number of filters', type=int)
parser.add_argument('--f0', required=True,
                    help='Central frequency of shortest filter', type=float)

parser.add_argument('--phases', nargs='+', default=['P', 'Pdiff'],
                    help='Phases to calculate kernels for')

parser.add_argument('--parameters', nargs='+', default=['vp'],
                    help='Model parameters to calculate kernels for')

parser.add_argument('--noplot', action="store_true", default=False,
                    help='Omit plotting (especially on headless environments)')

args = parser.parse_args()

model_params = args.parameters
nparam = len(model_params)

# Read event file
if args.event is None:
    time = obspy.core.UTCDateTime()
    origin = obspy.core.event.Origin(latitude=90.0, longitude=0.0,
                                     depth=00.0e3, time=time)

    focmec = obspy.core.event.FocalMechanism()
    MT = obspy.core.event.MomentTensor(
        tensor=obspy.core.event.Tensor(m_rr=0e20, m_pp=-1e20, m_tt=1e20,
                                       m_rp=0e20, m_rt=0e20, m_tp=1e20))
    focmec.moment_tensor = MT

    mag = obspy.core.event.Magnitude(mag=7.0)

    desc = obspy.core.event.EventDescription(text='Northpole Quake',
                                             type='earthquake name')

    event = obspy.core.event.Event(event_descriptions=[desc])
    event.magnitudes.append(mag)
    event.origins.append(origin)
    event.focal_mechanisms.append(focmec)

else:
    event = obspy.read_events(args.event)[0]

origin = event.origins[0]

event.write('CMTSOLUTION', format='CMTSOLUTION')

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
        fid.write('%s %s %4.1f 0.5 %4.1f 0\n' %
                  (filter_name[ifilter],
                   'Gabor',
                   filter_freq[ifilter],
                   filter_freq[ifilter]))

# Version to loop over stations in inventory
model = TauPyModel(model='ak135')

# Length of wavefield database
t_max = 1800.

if args.receivers is None:
    stats_sort = []
    dist_min = 20.0
    dist_max = 100.0
    dist_step = 20.0

    dists = np.arange(dist_min, dist_max + dist_step, dist_step)

    stats = []
    for irec in range(0, len(dists)):
        stat = obspy.core.inventory.Station(latitude=90.0 - dists[irec],
                                            longitude=00,
                                            code='R%03d' % int(dists[irec]),
                                            elevation=0)
        stats.append(stat)

else:
    # Read stations file
    inventory = obspy.read_inventory(args.receivers[0])
    for recfile in args.receivers[1:-1]:
        inventory += obspy.read_inventory(recfile)

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
stats_sort = sorted(stats, key=lambda stat:
                    locations2degrees(stat.latitude,
                                      stat.longitude,
                                      origin.latitude,
                                      origin.longitude))
nrec = len(stats_sort)

with open('receiver.dat', 'w') as fid:
    fid.write('%d   # nr of receivers\n' % nrec)
    fid.write('Z    # component\n')

    phase_list = args.phases

    for stat in stats_sort:
        distance = locations2degrees(stat.latitude, stat.longitude,
                                     origin.latitude, origin.longitude)
        # print distance
        tt = model.get_travel_times(distance_in_degree=distance,
                                    source_depth_in_km=origin.depth*1e-3,
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
                                                  nphase * args.nfilter *
                                                  nparam))
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
                                   travel_time - 10,
                                   travel_time + filter_length[iband]-10,
                                   param))

if not args.noplot:
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection=ccrs.Robinson())
    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()

    # ax.stock_img()
    ax.coastlines()

    ax.plot(origin.longitude, origin.latitude, '*', color='red',
            markersize=15, transform=ccrs.PlateCarree())

    # add a circle
    for dist in (30, 90, 120):
        circle = mpatches.Circle((origin.longitude, origin.latitude), dist,
                                 ec="black", fc="None", ls="dotted",
                                 transform=ccrs.PlateCarree())
        ax.add_patch(circle)

    # stats = []
    # stat_names = []
    # # Create list of stations, independent of network
    # for network in inventory:
    #     for stat in network:
    #         if (stat.code not in stat_names):
    #             stats.append(stat)
    #             stat_names.append(stat.code)
    #
    # # Sort stations by distance to source
    # stats_sort = sorted(stats[0:nrec], key=lambda
    #                     stat: locations2degrees(stat.latitude, stat.longitude,
    #                                             origin.latitude,
    #                                             origin.longitude))
    # nrec = len(stats_sort)

    for stat in stats_sort:
        ax.plot(stat.longitude, stat.latitude, 'v', color='darkgreen',
                transform=ccrs.PlateCarree(), markersize=10)
    ax.legend(fontsize=12, numpoints=1, markerscale=2.0)

    fname = './stations.png'
    fig.savefig(filename=fname,
                bbox_inches='tight',
                dpi=300)

# import numpy as np
# lat1 = stat.latitude * np.pi / 180
# lon1 = stat.longitude * np.pi / 180
# lat2 = event.latitude * np.pi / 180
# lon2 = event.longitude * np.pi / 180
# v1 = (np.cos(lon1)*np.cos(lat1), np.sin(lon1)*np.cos(lat1), np.cos(lat1))
# v2 = (np.cos(lon2)*np.cos(lat2), np.sin(lon2)*np.cos(lat2), np.cos(lat2))
#
# print np.cross(v1, v2)
