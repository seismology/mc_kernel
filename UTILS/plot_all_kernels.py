# coding: utf-8
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

from paraview import simple

import os
import numpy as np

# Create output directory
kernel_plot_dir = './Kernel_plots'
if not os.path.exists(kernel_plot_dir):
    os.mkdir(kernel_plot_dir)

# disable automatic camera reset on 'Show'
simple._DisableFirstRenderCameraReset()

# Reset render view, seems to change slightly in pvpython mode
renderView1 = simple.CreateView('RenderView')
renderView1.ViewSize = [960, 540]
renderView1.CameraPosition = [5103005.0, -34811145.0, 0.0]
renderView1.CameraFocalPoint = [5103005.0, 0.0, 0.0]
renderView1.CenterOfRotation = [253031.0, 0.0, -65653.0]
renderView1.CameraViewUp = [0.0, 0.0, 1.0]
renderView1.InteractionMode = '2D'
renderView1.OrientationAxesVisibility = 0
renderView1.OrientationAxesLabelColor = [0.313, 0.313, 0.313]
renderView1.StereoType = 0
renderView1.LightSwitch = 1
renderView1.LightIntensity = 0.2
renderView1.CameraParallelScale = 6520917.036707207
renderView1.Background = [1.0, 1.0, 1.0]


kerner_kernelxdmf = simple.XDMFReader(FileNames=['kerner_kernel.xdmf'])

kerner_kernelxdmf.GridStatus = ['grid']

kerner_all = simple.servermanager.Fetch(kerner_kernelxdmf)

# Apparently, it is possible to use any name here. We just need a
# transfer function which we can modify later
LUT = simple.GetColorTransferFunction('HalliGalli')
LUT.RGBPoints = [-1e-10, 0.231373, 0.298039, 0.752941,
                 0, 1.0, 1.0, 1.0,
                 1e-10, 0.705882, 0.0156863, 0.14902]
LUT.NanColor = [0.500008, 0.0, 0.0]
LUT.ScalarRangeInitialized = 1.0

# get any opacity transfer function/opacity map
# Not used so far
PWF = simple.GetOpacityTransferFunction('HalliGalli')
PWF.Points = [-1e-10, 0.0, 0.5, 0.0, 1e-10, 1.0, 0.5, 0.0]
PWF.ScalarRangeInitialized = 1

# setup the color legend parameters for each legend in this view

# get color legend/bar for LUT in view renderView1
LUTColorBar = simple.GetScalarBar(LUT, renderView1)
LUTColorBar.Position = [0.8322310304209063, 0.50956937799043066]
LUTColorBar.Position2 = [0.12, 0.43000000000000077]
LUTColorBar.Title = 'K_x [s/m^3]'
LUTColorBar.ComponentTitle = ''
LUTColorBar.TitleColor = [0.3137, 0.3137, 0.3137]
LUTColorBar.LabelColor = [0.3137, 0.3137, 0.3137]

# Loop over all Point variables in file and print them
for data in kerner_kernelxdmf.PointData:
    if data.Name[0:3] == 'K_x':

        # Get complete array of this variable
        npoints = kerner_all.GetNumberOfPoints()
        data_list = []
        for i in range(0, npoints):
            newval = kerner_all.GetPointData().GetArray(data.Name).GetValue(i)
            data_list.append(newval)
        data_array = np.asarray(data_list)

        # Set color range to 99 percentile
        maxval = np.percentile(abs(data_array), 99.5)
        maxval = round(maxval, int(-np.log10(maxval))+1)

        LUT.RGBPoints = [-maxval, 0.231373, 0.298039, 0.752941,
                         0, 1.0, 1.0, 1.0,
                         maxval, 0.705882, 0.0156863, 0.14902]

        print 'Plotting kernel %s, range(%f, %f)' % (data.Name,
                                                     maxval, maxval)

        # show data from kerner_kernelxdmf
        kerner_kernelxdmfDisplay = simple.Show(kerner_kernelxdmf, renderView1)

        # trace defaults for the display properties.
        kerner_kernelxdmfDisplay.ColorArrayName = ['POINTS', data.Name]
        kerner_kernelxdmfDisplay.LookupTable = LUT

        # show color legend
        kerner_kernelxdmfDisplay.SetScalarBarVisibility(renderView1, True)

        filename_out = os.path.join(kernel_plot_dir, '%s.png' % data.Name)
        simple.SaveScreenshot(filename=filename_out,
                              view=renderView1,
                              magnification=2)

        print '  ...done!'

# Loop over all Cell variables in file and print them
for data in kerner_kernelxdmf.CellData:
    if data.Name[0:3] == 'K_x':

        # Get complete array of this variable
        ncells = kerner_all.GetNumberOfCells()
        data_list = []
        for i in range(0, ncells):
            newval = kerner_all.GetCellData().GetArray(data.Name).GetValue(i)
            data_list.append(newval)
        data_array = np.asarray(data_list)

        # Set color range 99 percentile
        maxval = np.percentile(abs(data_array), 99.5)
        maxval = round(maxval, int(-np.log10(maxval))+1)

        LUT.RGBPoints = [-maxval, 0.231373, 0.298039, 0.752941,
                         0, 1.0, 1.0, 1.0,
                         maxval, 0.705882, 0.0156863, 0.14902]

        print 'Plotting kernel %s, range(%f, %f)' % (data.Name,
                                                     maxval, maxval)

        # show data from kerner_kernelxdmf
        kerner_kernelxdmfDisplay = simple.Show(kerner_kernelxdmf, renderView1)

        # trace defaults for the display properties.
        kerner_kernelxdmfDisplay.ColorArrayName = ['CELLS', data.Name]
        kerner_kernelxdmfDisplay.LookupTable = LUT

        # show color legend
        kerner_kernelxdmfDisplay.SetScalarBarVisibility(renderView1, True)

        filename_out = os.path.join(kernel_plot_dir, '%s.png' % data.Name)
        simple.SaveScreenshot(filename=filename_out,
                              view=renderView1,
                              magnification=2)

        print '  ...done!'
