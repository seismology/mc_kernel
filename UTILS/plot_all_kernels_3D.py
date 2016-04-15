# coding: utf-8
#******************************************************************************
#
#    This file is part of:
#    MC Kernel: Calculating seismic sensitivity kernels on unstructured meshes
#    Copyright (C) 2016 Simon Staehler, Martin van Driel, Ludwig Auer
#
#    You can find the latest version of the software at:
#    <https://www.github.com/tomography/mckernel>
#
#    MC Kernel is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    MC Kernel is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with MC Kernel. If not, see <http://www.gnu.org/licenses/>.
#
#******************************************************************************

"""
Plot all kernels with Paraview
Author: Simon Stähler, LMU München

This script uses Paraview, so it should be started using pvpython (>=4.3)
The kernel pictures can later be combined with seismogram pictures created
with plot_all_seismograms.py, using the script create_composite_plots.py

The camera focal point is adapted to the 2D Mesh in Meshes/circle_050.inp

To plot 3D kernels, you have to modify this script by introducing a clipping
plane at the right place.
"""

from paraview.simple import *
import numpy as np
import os

def read_receiver_dat(src_file, rec_file):

    with open(src_file) as f:
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        event_lat = float(f.readline().split()[1])
        event_lat *= np.pi / 180
        event_lon = float(f.readline().split()[1])
        event_lon *= np.pi / 180

    with open(rec_file) as f:
        nkernel_total = 0
        fullstrain_kernel = False
        # Read number of receivers
        str_line = f.readline()
        nrec = int(str_line.split()[0])
        print 'Number of receivers: %d' % nrec
        # Read seismogram component
        str_line = f.readline()
        # seis_cmp = str_line.split()[0]
        for irec in range(0, nrec):
            str_line = f.readline()
            rec_name = str_line.split()[0]
            rec_lat = float(str_line.split()[1])
            rec_lon = float(str_line.split()[2])
            nkernel = int(str_line.split()[3])
            print 'Receiver: %s, coordinates: (%f, %f), %d kernels' % \
                  (rec_name, rec_lat, rec_lon, nkernel)
            for ikernel in range(0, nkernel):
                str_line = f.readline()
                kernel_name = str_line.split()[0]
                # filter_name = str_line.split()[1]
                # misfit_name = str_line.split()[2]
                #time_window_start = float(str_line.split()[3])
                #time_window_stop = float(str_line.split()[4])

                #kernel_var_name ='K_x_%s_%s_%03d' % (rec_name, kernel_name,
                #                                     ikernel)
                kernel_var_name ='K_x_%s_%s' % (rec_name, kernel_name)

                lat1 = rec_lat * np.pi / 180
                lon1 = rec_lon * np.pi / 180
                v1 = (np.cos(lon1)*np.cos(lat1),
                      np.sin(lon1)*np.cos(lat1),
                      np.sin(lat1))
                v2 = (np.cos(event_lon)*np.cos(event_lat),
                      np.sin(event_lon)*np.cos(event_lat),
                      np.sin(event_lat))

                normal_vector = np.cross(v1, v2)

                # Normalize normal vector
                normal_vector /= np.sqrt(np.sum(normal_vector**2))

                try:
                    print_point(kernel_var_name, normal_vector)
                except:
                    print_cell(kernel_var_name, normal_vector)

            nkernel_total += nkernel

        print 'Number of kernels: %d' % nkernel_total

def print_point(name, normalvector):

    renderView1.CameraPosition = normalvector * 3e7 #+ [2e7, 0, 0]
    renderView1.CameraFocalPoint = [0.0, 0.0, 0.0]

    # create a new 'Clip'
    clip1 = Clip(Input=kerner_kernelxdmf)
    clip1.ClipType = 'Plane'
    clip1.Scalars = ['POINTS', name]
    clip1.Value = -3.23037052154541
    clip1.InsideOut = 1
    clip1.Crinkleclip = 1

    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [0.0, 0.0, 0.0]
    clip1.ClipType.Normal = normalvector

    # Set color range to 10% of maximum value 
    # Get maxval - not elegant
    for data in clip1.PointData:
        if data.Name == name:
            maxval = np.max(-data.GetRange()[0], data.GetRange()[1]) * 0.1
            maxval = round(maxval, int(-np.log10(maxval))+1)

            LUT.RGBPoints = [-maxval, 0.231373, 0.298039, 0.752941,
                             0, 1.0, 1.0, 1.0,
                             maxval, 0.705882, 0.0156863, 0.14902]

    print 'Plotting kernel %s, range(%f, %f)' % (name, maxval, maxval)


    # show data from kerner_kernelxdmf
    ClipDisplay = Show(clip1, renderView1)

    # trace defaults for the display properties.
    ClipDisplay.ColorArrayName = ['POINTS', name]
    ClipDisplay.LookupTable = LUT

    # show color legend
    ClipDisplay.SetScalarBarVisibility(renderView1, True)

    filename_out = os.path.join(kernel_plot_dir, '%s.png' % name)
    SaveScreenshot(filename=filename_out,
                   view=renderView1,
                   magnification=2)

    print '  ...done!'
    Delete(clip1)

# Loop over all Cell variables in file and print them
def print_cell(name, normalvector):

    renderView1.CameraPosition = normalvector * 3e7 #+ [2e7, 0, 0]
    renderView1.CameraFocalPoint = [0.0, 0.0, 0.0]

    # create a new 'Clip'
    clip1 = Clip(Input=kerner_kernelxdmf)
    clip1.ClipType = 'Plane'
    clip1.Scalars = ['CELLS', name]
    clip1.Value = -3.23037052154541
    clip1.InsideOut = 1
    clip1.Crinkleclip = 1

    # Set color range to 10% of maximum value 
    # Get maxval - not elegant
    for data in clip1.CellData:
        if data.Name == name:
          maxval = np.max(np.abs(data.GetRange())) * 0.1
          maxval = round(maxval, int(-np.log10(maxval))+1)

          LUT.RGBPoints = [-maxval, 0.231373, 0.298039, 0.752941,
                           0, 1.0, 1.0, 1.0,
                           maxval, 0.705882, 0.0156863, 0.14902]

    print 'Plotting kernel %s, range(%f, %f)' % (name, maxval, maxval)


    # init the 'Plane' selected for 'ClipType'
    clip1.ClipType.Origin = [0.0, 0.0, 0.0]
    clip1.ClipType.Normal = normalvector

    # show data from kerner_kernelxdmf
    ClipDisplay = Show(clip1, renderView1)

    # trace defaults for the display properties.
    ClipDisplay.ColorArrayName = ['CELLS', name]
    ClipDisplay.LookupTable = LUT

    # show color legend
    ClipDisplay.SetScalarBarVisibility(renderView1, True)

    filename_out = os.path.join(kernel_plot_dir, '%s.png' % name)
    SaveScreenshot(filename=filename_out,
                   view=renderView1,
                   magnification=2)

    print '  ...done!'

    Delete(clip1)



# Create output directory
kernel_plot_dir = './Kernel_plots'
if not os.path.exists(kernel_plot_dir):
    os.mkdir(kernel_plot_dir)


# disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Reset render view, seems to change slightly in pvpython mode
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [960, 540]
renderView1.CameraPosition = [5103005.0, -34811145.0, 0.0]
renderView1.CameraFocalPoint = [5103005.0, 0.0, 0.0]
renderView1.CenterOfRotation = [0.0, 0.0, 0.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
#renderView1.InteractionMode = '2D'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.313, 0.313, 0.313]
renderView1.StereoType = 0
renderView1.LightSwitch = 1
renderView1.LightIntensity = 0.2
renderView1.CameraParallelScale = 6520917.036707207
renderView1.Background = [1.0, 1.0, 1.0]

# create a new 'Arrow'
arrow1 = Arrow()
arrow1.TipResolution = 64
arrow1.TipRadius = 0.01
arrow1.TipLength = 0.025
arrow1.ShaftResolution = 64
arrow1.ShaftRadius = 0.005

# show data from arrow1
arrow1Display = Show(arrow1, renderView1)

# trace defaults for the display properties.
arrow1Display.ColorArrayName = [None, '']
arrow1Display.DiffuseColor = [1.0, 0.0, 0.0]
arrow1Display.Position = [0.0, 0.0, -10000000.0]
arrow1Display.Scale = [20000000.0, 20000000.0, 20000000.0]
arrow1Display.Orientation = [0.0, 270.0, 0.0]

# create a new 'Sphere'
# Core
sphere1 = Sphere()
sphere1.Radius = 3461000.0
sphere1.ThetaResolution = 256
sphere1.PhiResolution = 256

# show data from sphere1
sphere1Display = Show(sphere1, renderView1)
# trace defaults for the display properties.
sphere1Display.ColorArrayName = [None, '']
sphere1Display.DiffuseColor = [0.4, 0.4, 0.4]

kerner_kernelxdmf = XDMFReader(FileNames=['kerner_kernel.xdmf'])

kerner_kernelxdmf.GridStatus = ['grid']

# Apparently, it is possible to use any name here. We just need a
# transfer function which we can modify later
LUT = GetColorTransferFunction('HalliGalli')
LUT.RGBPoints = [-1e-10, 0.231373, 0.298039, 0.752941,
                 0, 1.0, 1.0, 1.0,
                 1e-10, 0.705882, 0.0156863, 0.14902]
LUT.LockScalarRange = 1
LUT.NanColor = [0.500008, 0.0, 0.0]
LUT.ScalarRangeInitialized = 1.0

# get any opacity transfer function/opacity map
# Not used so far
PWF = GetOpacityTransferFunction('HalliGalli')
PWF.Points = [-1e-10, 0.0, 0.5, 0.0, 1e-10, 1.0, 0.5, 0.0]
PWF.ScalarRangeInitialized = 1

# setup the color legend parameters for each legend in this view

# get color legend/bar for LUT in view renderView1
LUTColorBar = GetScalarBar(LUT, renderView1)
LUTColorBar.Position = [0.8322310304209063, 0.50956937799043066]
LUTColorBar.Position2 = [0.12, 0.43000000000000077]
LUTColorBar.Title = 'K_x [s/m^3]'
LUTColorBar.ComponentTitle = ''
LUTColorBar.TitleColor = [0.3137, 0.3137, 0.3137]
LUTColorBar.LabelColor = [0.3137, 0.3137, 0.3137]

# Load Coastlines
coastlinevtk = LegacyVTKReader(FileNames=['../../UTILS/Coastlines/Coastline.vtk'])
# show data from coastlinevtk
coastlinevtkDisplay = Show(coastlinevtk, renderView1)
# trace defaults for the display properties.
coastlinevtkDisplay.ColorArrayName = ['POINTS', '']
coastlinevtkDisplay.DiffuseColor = [0.2, 0.2, 0.2]
coastlinevtkDisplay.Scale = [1000.0, 1000.0, 1000.0]

read_receiver_dat('CMTSOLUTION', 'receiver.dat')

