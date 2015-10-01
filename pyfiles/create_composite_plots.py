"""
Create a composite plot of kernel and seismogram
Author: Simon Stähler, LMU München

Before running this script, it is necessary to create kernel pictures with
plot_all_kernels.py and seismogram pictures with plot_all_seismograms.py
"""
# coding: utf-8
from wand.image import Image
import os
import glob

seis_plot_dir = './Seismogram_plots'
kernel_plot_dir = './Kernel_plots/'
comp_plot_dir = './Composite_plots/'

if not os.path.exists(comp_plot_dir):
    os.mkdir(comp_plot_dir)

# Get list of all seismogram plot files
smgr_list = glob.glob(os.path.join(seis_plot_dir, '*.png'))

for filename_smgr in smgr_list:
    kernel = os.path.split(filename_smgr)[1]
    # Just assume that there is a plot file for each kernel
    filename_kernel = os.path.join(kernel_plot_dir, 'K_x_%s' % kernel)

    filename_comp = os.path.join(comp_plot_dir, 'Composite_%s.png' % kernel)
    print 'Processing Kernel: %s...' % kernel

    # Load kernel plot and seismogram plot
    im_kernel = Image(filename=filename_kernel)
    im_smgr = Image(filename=filename_smgr)
    im_smgr.crop(left=100)

    # Create composite plot
    im_kernel.composite(im_smgr, left=1100, top=580)

    im_kernel.save(filename=filename_comp)
    print '  ...done!'
