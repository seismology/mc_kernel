import argparse
import os
import shutil
import glob
import datetime
import subprocess
from netCDF4 import Dataset 

"""
Find 2^n that is equal to or greater than.
"""
def nextpow2(i):
  n = 1
  while n < i: n *= 2
  return n

def estimate_memory():
  # Estimate memory usage which cannot be controlled by input parameters

  memory_mesh = 4 * ( 2*npoints_fwd + 4*nelems_fwd + 4*nelems_fwd + 25*nelems_fwd) * 2
  print 'Mesh size in memory: %f MB'%(memory_mesh/(2**20))

  nomega = nextpow2(ndumps_fwd)
  if full_strain:
    ndim = 6
  else:
    ndim = 1

  memory_fft = 8 * (nomega + ndumps_fwd) * float(params['ELEMENTS_PER_TASK']) * ndim
  print 'Memory for FFT types: %f MB'%(memory_fft/(2**20))

  memory_fields = 8 * float(params['ELEMENTS_PER_TASK']) * ndim * (5*nomega + 2*ndumps_fwd)
  print 'Memory for Wavefields: %f MB'%(memory_fields/(2**20))

  # The following memory requirements were determined using the Massif heap profiler
  # from the Valgrind package on two Ubuntu 14.04 machines. It might be that the
  # memory requirements of the HDF5 library are different on other architectures
  # Memory usage of HDF5 library seems to be pretty constant around 120 MB.
  memory_hdf5 = 150 * 2**20 
  print 'Memory for HDF5/NetCDF4 library: %f MB'%(memory_hdf5/(2**20))

  # Memory usage of KD-Trees is roughly 80 Byte per mesh point
  memory_kdtree = 80 * npoints_fwd  
  print 'Memory for KD-Trees: %f MB'%(memory_kdtree/(2**20))

  return memory_mesh + memory_fft + memory_hdf5 + memory_kdtree + memory_fields



def auto_buffer_size(memory_available):
  if full_strain:
    ndim = 6
  else:
    ndim = 1

  memory_for_buffers = (memory_available - estimate_memory())*0.9

  # Rule: Strain buffer gets 90% of the available memory, displ. buffer 10%
  size_strain_buffer = int(memory_for_buffers * 0.9 / (25 * 4 * ndumps_fwd * 6 * ndim))
  size_disp_buffer   = int(memory_for_buffers * 0.1 / (3 * 4 * ndumps_fwd * 6 ))

  memory_buffers_strain = 25 * 4 * ndumps_fwd * 6 * size_strain_buffer * ndim
  print 'Strain buffer size: %f MB'%(memory_buffers_strain/(2**20))
  memory_buffers_disp = 3 * 4 * ndumps_fwd * 6 * size_disp_buffer
  print 'Displ. buffer size: %f MB'%(memory_buffers_disp/(2**20))

  
  return str(size_strain_buffer), str(size_disp_buffer)

def read_receiver_dat(rec_file):
  with open(rec_file) as f:

    nkernel_total = 0
    fullstrain_kernel = False
    # Read number of receivers
    str_line = f.readline()
    nrec = int(str_line.split()[0])
    print 'Number of receivers: %d'%nrec

    # Read seismogram component
    str_line = f.readline()
    seis_cmp = str_line.split()[0]
    #print 'Seismogram component: %s'%seis_cmp

    for irec in range(0, nrec):
        str_line = f.readline()
        rec_name = str_line.split()[0]
        rec_lat  = float(str_line.split()[1])
        rec_lon  = float(str_line.split()[2])
        nkernel  = int(str_line.split()[3])
        
        #print 'Receiver: %s, coordinates: (%f, %f), %d kernels'%(rec_name, rec_lat, rec_lon, nkernel)
        
        for ikernel in range(0, nkernel):
            str_line = f.readline()
            #print str_line
            kernel_name = str_line.split()[0]
            filter_name = str_line.split()[1]
            misfit_name = str_line.split()[2]
            time_window_start = float(str_line.split()[3])
            time_window_stop  = float(str_line.split()[4])
            model_param = str_line.split()[5]
            if model_param in ('vs', 'rho', 'vsh', 'vsv', 'eta', 'phi', 'xi', 'mu '):
              fullstrain_kernel = True
            elif model_param not in ('lam', 'vp', 'vph', 'vpv'):
              raise RuntimeError('Unknown model parameter %s in %s'%(model_param, rec_file))

        nkernel_total += nkernel
    
    print 'Number of kernels: %d'%nkernel_total

  return nkernel_total, nrec, fullstrain_kernel

def define_arguments():
  parser = argparse.ArgumentParser(description='Create Kerner input file and submit job.',
                                   formatter_class=argparse.RawTextHelpFormatter)

  parser.add_argument('job_name', help='Job directory name')

  helptext = """Input file to use. It will overwrite default values, 
  but will be overwritten by any argument to this function."""
  parser.add_argument('-i', '--input_file', help=helptext)

  parser.add_argument('-n', '--nslaves', type=int, default=2,
                      metavar='N', 
                      help='Number of slaves to use')

  parser.add_argument('-m', '--message', metavar='JOB_DESCRIPTION_MESSAGE',
                      help="Description of run, which is saved in job_name/README.run\n"+
                           "If omitted, an editor window opens to collect description.")

  parser.add_argument('--what_to_do', choices=['integrate_kernel', 'plot_wavefield'], 
                      default='integratekernel',
                      help='Calculate kernels or just plot wavefields')

  parser.add_argument('-a', '--available_memory', type=int,
                      help='Amount of memory available in MB')

  parser.add_argument('-q', '--queue', choices=['SuperMUC', 'local', 'monch'],
                      default='local',
                      help='Queue to use. Default is local, which starts a job with MPIRUN')

  ###############################################################################
  # Queue options
  ###############################################################################
  hpc_queue = parser.add_argument_group('Options specific to HPC queues')
  hpc_queue.add_argument('--wall_time', type=int, metavar='WALLTIME_IN_H',
                         default=10,
                         help='Walltime in hours')
  hpc_queue.add_argument('--mail_address', help='Mail address for HPC notifications')


  ###############################################################################
  # AxiSEM wavefield directories
  ###############################################################################

  axisem_dirs = parser.add_argument_group('AxiSEM run directories')
  axisem_dirs.add_argument('--fwd_dir', default='./wavefield/fwd/',
                      help='Path to AxiSEM forward run')
  axisem_dirs.add_argument('--bwd_dir', default='./wavefield/bwd/',
                      help='Path to AxiSEM backward run')

  ###############################################################################
  # input files
  ###############################################################################

  input_files = parser.add_argument_group('Required input files')
  input_files.add_argument('--src_file', default='CMTSOLUTION',
                      help='Path to source solution file in CMTSOLUTION format')
  input_files.add_argument('--rec_file', default='receiver.dat',
                      help='Path to receiver and kernel file')
  input_files.add_argument('--filt_file', default='filters.dat',
                      help='Path to filter file')
  input_files.add_argument('--stf_file', default='stf_20s.dat',
                      help='Path to Source Time Function file')

  ###############################################################################
  # Mesh file-related options
  ###############################################################################

  mesh_files = parser.add_argument_group('Inversion mesh')
  mesh_helptext = """
  Select the mesh file type. Allowed values are 
  abaqus      : .inp file, can be generated with Qubit or other codes. Can
                contain various geometries and multiple sub-objects
                Supported geometries: tetrahedra, triangles, quadrilaterals
                Set file name in MESH_FILE_ABAQUS 

  tetrahedral : tetrahedral mesh in two separate files with 
                1. coordinates of the vertices (MESH_FILE_VERTICES)
                2. the connectivity of the facets of the tetrahedrons
                   (MESH_FILE_FACETS)"""
  mesh_files.add_argument('--mesh_file_type', default='abaqus',
                          choices=['abaqus', 'tetrahedral'], help=mesh_helptext)
  mesh_files.add_argument('--mesh_file_abaqus', default='Meshes/flat_triangles.inp',
                      help='Path to Abaqus mesh file')
  mesh_files.add_argument('--mesh_file_vertices', default='unit_tests/vertices.TEST',
                      help='Path to Vertices file (only if --mesh_file_type=tetrahedral)')
  mesh_files.add_argument('--mesh_file_facets', default='unit_tests/facets.TEST', 
                      help='Path to Facets file (only if --mesh_file_type=tetrahedral)')

  ###############################################################################
  # Kernel-related options
  ###############################################################################

  kernel_options = parser.add_argument_group('Kernel calculation options')
  helptext = """
  Calculate kernels for absolute values Vp instead of relative perturbations dVp 
  with respect to the background model""" 

  kernel_options.add_argument('--kernel_for_absolute_perturbations', action="store_true", default=False,
                              help=helptext)

  helptext = """On which base functions are the kernels defined?
  volumetric (default): Each voxel is a base function (Boschi & Auer)
  onvertices:           Each vertex has a set of non-orthogonal base functions
                        defined on it (Nolet & Sigloch)"""
  kernel_options.add_argument('--int_type', choices=['volumetric', 'onvertices'], 
                              default='volumetric', help=helptext)
  helptext = """
  For plotting reasons one may wish to skip the integration over cell-volume.
  Resulting kernels bear the unit [s/m^3]"""
  kernel_options.add_argument('--no_int_over_volume', action="store_true", default=False,
                              help=helptext)


  ###############################################################################
  # Monte Carlo-related options
  ###############################################################################

  mc_options = parser.add_argument_group('Monte Carlo options')
  helptext = """
  Number of points on which the kernel is evaluated per 
  Monte Carlo step. Default value is 4."""
  mc_options.add_argument('--points_per_mc_step', type=int, default=4,
                          help=helptext)
  helptext = """Maximum number of Monte Carlo iterations. Allows to skip 
  evaluation of single problematic cells. Default value is 1E6"""
  mc_options.add_argument('--maximum_iterations', type=int, default=1000000,
                          help=helptext)
  helptext = """Allowed absolute error before Monte Carlo integration is considered
  to be converged. When calculating this value, the volume is not considered,
  no matter whether --no_int_over_volume is set or not."""
  mc_options.add_argument('--allowed_error', type=float, default=1e-4,
                          help=helptext)
  helptext = """Allowed relative error before Monte Carlo integration in one cell
  is considered to be converged. Default value is 1e-2"""
  mc_options.add_argument('--allowed_relative_error', type=float, default=1e-2,
                          help=helptext)
  helptext = """Use pseudorandom numbers instead of quasirandom"""
  mc_options.add_argument('--use_pseudorandom_numbers', action="store_true", default=False,
                          help=helptext)

  ###############################################################################
  # Debugging-related options
  ###############################################################################

  debug_options = parser.add_argument_group('Debugging options')
  helptext = """
  This activates the optional (linearity test) to integrate relative kernels over
  model perturbations in percent, to assess how well our kernels predict measured
  traveltime perturbations for the same model. This only makes sense when not 
  calculating absolute kernels. """
  debug_options.add_argument('--int_over_3d_heterogeneities', action="store_true", default=False,
                      help=helptext)
  helptext = """Path to heterogeneity file"""
  debug_options.add_argument('--het_file', default='unit_tests/savani.rtpv',
                      help=helptext)

  helptext = """
  Integrate the kernel over the background model. Classically, this was assumed
  to result in the travel time of a phase. This assumption is highly dubious for
  wavefield-derived kernels. For legacy reasons, we can still leave it in. 
  Adds a version of the background model interpolated on the inversion mesh to 
  the output file."""
  debug_options.add_argument('--int_over_background_model', action="store_true", default=False,
                      help=helptext)

  helptext = """ Every slave writes out the values of all the kernels and their respective 
  estimated errors into his OUTPUT_??? file after each MC step. This can lead 
  to huge ASCII files (>1GB) with inane line lengths (approx. 20 x nkernel).
  However, it might be interesting to study the convergence behaviour. """
  debug_options.add_argument('--write_detailed_convergence', action="store_true", default=False,
                             help=helptext)

  helptext = """Do not deconvolve the Source Time Function and reconvolve with 
  the one set in --stf_file, but just timeshift the wavefields."""
  debug_options.add_argument('--no_deconvolve_stf', action="store_true", default=False,
                             help=helptext)

  helptext = """Integration scheme to calculate scalar kernels from seismograms and waveforms. 
  parseval (default):  Integration in frequency domain, using the Parseval theorem.
  trapezoidal:         Integration in time domain using the trapezoidal rule."""
  debug_options.add_argument('--integration_scheme', choices=['parseval', 'trapezoidal'], 
                             default='parseval', help=helptext)


  ###############################################################################
  # Output-related options
  ###############################################################################

  output_options = parser.add_argument_group('Output options')
  helptext = """Output format when dumping kernels and wavefields. 
  Choose between xdmf, Yale-style csr binary format (compressed sparse row) and
  ascii. Yet, the allowed error below is assumed as the truncation threshold in 
  csr and ascii storage"""
  output_options.add_argument('--dump_type', choices=['xdmf', 'ascii', 'csr'], default='xdmf',
                      help=helptext)
  helptext = """Write out Seismograms (raw full trace, filtered full trace 
  and cut trace) into run_dir/SEISMOGRAMS. Produces three files per kernel. 
  Disable to avoid congesting your file system."""
  output_options.add_argument('--write_seismograms', default=False,
                              help=helptext)
  helptext = """Prefix of output file names.
  Kernel files are called $OUTPUT_FILE_kernel.xdmf
  Wavefield movies are called $OUTPUT_FILE_wavefield.xdmf"""
  output_options.add_argument('--out_prefix', default='kerner',
                              help=helptext)

  ###############################################################################
  # Performance-related options
  ###############################################################################

  performance_options = parser.add_argument_group('Performance-related options')
  helptext = """Size of buffer for strain values. Since the strain has to be 
  calculated from the displacement stored in the AxiSEM files,
  increasing this buffer size saves IO access and CPU time."""
  performance_options.add_argument('--strain_buffer_size', type=int, default=1000,
                              help=helptext)

  helptext = """Size of buffer for displacement values. Displacement values are
  use to calculate strain later. Having a separate buffer here allows to save 
  some IO accesses."""
  performance_options.add_argument('--displ_buffer_size', type=int, default=100,
                              help=helptext)
  helptext = """Number of elements per Slave task. A larger value allows to the
  Slave to have more contiguous parts of the earth to work on, smaller values
  improve load balancing. It should be chosen such that each slave gets at least
  50-100 tasks to work on."""
  performance_options.add_argument('--elements_per_task', type=int, default=100,
                              help=helptext)
  helptext = """Do not sort the mesh elements. Just for debugging."""
  performance_options.add_argument('--no_sort_mesh_elements', action="store_true", default=False,
                              help=helptext)
  helptext = """Mask the source and the receiver element and set the kernel to 
  zero in each. A rough way to avoid spending hours until convergence in these 
  two elements in reached."""
  performance_options.add_argument('--mask_source_receiver', action="store_true", default=False,
                              help=helptext)
  helptext = """Dampen the kernel in a radius around source and receiver. 
  If a negative value is chosen, damping is switched off (default)."""
  performance_options.add_argument('--damp_radius_source_receiver', type=float, default=-100.E3,
                                   help=helptext)
  helptext = """FFTW Planning to use 
  Options: 
  ESTIMATE:   Use heuristic to find best FFT plan
  MEASURE:    Compute several test FFTs to find best plan (default)
  PATIENT:    Compute a lot of test FFTs to find best plan
  EXHAUSTIVE: Compute an awful amount of test FFTs to find best plan
  This option did not prove to be very useful on most systems."""
  performance_options.add_argument('--fftw_plan', default='MEASURE',
                      choices=['ESTIMATE', 'MEASURE', 'PATIENT', 'EXHAUSTIVE'],
                      help=helptext)
  return parser

parser = define_arguments()

args = parser.parse_args()

# Parse input file, if one was given
if args.input_file: 
  with open(args.input_file) as f:
    args_input_file = {}
    for line in f:
      # Skip comment and empty lines
      if line[0]!='#' and line.strip()!='':
        (key, val) = line.split()
        args_input_file[key] = val

# Merge input variables from input_file, arguments and default values
params = {}
# Loop over all possible arguments
for key, value in vars(args).iteritems():
  if not key in ('nslaves', 'job_name', 'queue', 'available_memory'): 
    # If an input file is selected, get values from there by default
    if args.input_file:
      if value == parser.get_default(key):
        key_out   = key.upper()
        try:
          value_out = str(args_input_file[key.upper()]).strip("'")
        except KeyError:
          # If value is not set in input file
          value_out = str(value)
      else:
        # Unless values were explicitly given
        key_out = key.upper()
        value_out = str(value)
    else:
      # In all other cases, take default values
      key_out = key.upper()
      value_out = str(value)

  params[key_out] = value_out

# Read receiver file and get number of receivers, kernels and whether the full strain
# has to be read for any kernel (increases the memory footprint of the buffers by factor of 6)
nrec, nkernel, full_strain = read_receiver_dat(args.rec_file)

# Check for AxiSEM wavefield files to get mesh size
fwd_path = os.path.join(os.path.realpath(params['FWD_DIR']), 'MZZ', 'Data', 'ordered_output.nc4')
bwd_path = os.path.join(os.path.realpath(params['BWD_DIR']), 'PZ', 'Data', 'ordered_output.nc4')
if os.path.exists(fwd_path):
  nc_fwd = Dataset(fwd_path)
  npoints_fwd = getattr(nc_fwd, "npoints")
  nelems_fwd  = getattr(nc_fwd, "nelem_kwf_global")
  ndumps_fwd  = getattr(nc_fwd, "number of strain dumps")
  nc_fwd.close()
else: 
  errmsg = 'Could not find a wavefield file in the fwd_dir %s\n%s'%(params['FWD_DIR'], fwd_path)
  raise IOError(errmsg)
  
if os.path.exists(bwd_path):
  nc_bwd = Dataset(bwd_path)
  npoints_bwd = getattr(nc_bwd, "npoints")
  nelems_bwd  = getattr(nc_bwd, "nelem_kwf_global")
  ndumps_bwd  = getattr(nc_bwd, "number of strain dumps")
  nc_bwd.close()
else: 
  errmsg = 'Could not find a wavefield file in the bwd_dir %s'%params['FWD_DIR']
  raise IOError(errmsg)


# Define buffer sizes based on available memory
if args.available_memory:
  params['STRAIN_BUFFER_SIZE'], params['DISPL_BUFFER_SIZE'] = auto_buffer_size(args.available_memory*2**20)
  print params['STRAIN_BUFFER_SIZE'], params['DISPL_BUFFER_SIZE'] 
  

# Sanity check, whether fwd and bwd mesh have the same sizes and the same number 
# of wavefield time steps.
if npoints_fwd!=npoints_bwd or nelems_fwd!=nelems_bwd or ndumps_fwd!=ndumps_bwd:
  raise RuntimeError('Forward and backward run did not use the same parameters')


params_out = {}

# Create run_dir
run_dir = args.job_name
os.mkdir(run_dir)

# Copy necessary files to rundir
for key, value in params.iteritems():

  if key in ('NSLAVES', 'JOB_NAME', 'MESH_FILE_ABAQUS', 'MESH_FILE_VERTICES', 'MESH_FILE_FACETS', 
             'HET_FILE', 'INPUT_FILE', 'MESSAGE', 'AVAILABLE_MEMORY', 'WALL_TIME', 'MAIL_ADDRESS'): 
    # Effectively nothing to do
    print ''
  elif key=='SRC_FILE':
    shutil.copy(value, os.path.join(run_dir, 'CMTSOLUTION'))
    params_out[key] = 'CMTSOLUTION'
  elif key=='REC_FILE':
    shutil.copy(value, os.path.join(run_dir, 'receiver.dat'))
    params_out[key] = 'receiver.dat'
  elif key=='FILT_FILE':
    shutil.copy(value, os.path.join(run_dir, 'filters.dat'))
    params_out[key] = 'filters.dat'
  elif key=='STF_FILE':
    shutil.copy(value, os.path.join(run_dir, 'stf.dat'))
    params_out[key] = 'stf.dat'

  elif key=='MESH_FILE_TYPE':
    params_out[key] = value
    if value=='abaqus':
      shutil.copy(params['MESH_FILE_ABAQUS'], os.path.join(run_dir, 'mesh.inp'))
      params_out['MESH_FILE_ABAQUS'] = 'mesh.inp'
    elif value=='tetrahedral':
      shutil.copy(params['MESH_FILE_VERTICES'], os.path.join(run_dir, 'mesh.VERTICES'))
      shutil.copy(params['MESH_FILE_FACETS'], os.path.join(run_dir, 'mesh.FACETS'))
      params_out['MESH_FILE_VERTICES'] = 'mesh.VERTICES'
      params_out['MESH_FILE_FACETS'] = 'mesh.FACETS'
    else:
      raise KeyError('Unknown mesh file type: %s'%args.mesh_file_type)

  elif key=='INT_OVER_3D_HETEROGENEITIES':
    params_out[key] = value
    shutil.copy(params['HET_FILE'], os.path.join(run_dir, 'heterogeneities.dat'))
    params_out['HET_FILE'] = 'heterogeneities.dat'
  elif key in ('FWD_DIR', 'BWD_DIR'):
    # Set mesh dir to absolute path
    params_out[key] = os.path.realpath(value)

  else:
    params_out[key] = value

# Open editor window to write run descriptor
out_readme = 'readme_temp.txt' 
f_readme = open(out_readme, 'w')
f_readme.write('MC KERNEL run for %d CPUs, started on %s\n'%(args.nslaves, 
                                                           str(datetime.datetime.now())))
f_readme.write('  by user ''%s'' on ''%s''\n'%(os.environ.get('USER'), 
                                               os.environ.get('HOSTNAME')))
if args.message:
  f_readme.write(args.message)
  f_readme.close()
else:
  f_readme.close()
  editor = os.environ.get('EDITOR')
  subprocess.check_call('%s %s'%(editor, out_readme), shell=True)

# Move README file to rundir
shutil.move(out_readme, os.path.join(run_dir, 'README.run'))

# Create directory for seismogram output
os.mkdir(os.path.join(run_dir, 'Seismograms'))

# Create input file for run
out_input_file = os.path.join(run_dir, 'inparam')
with open(out_input_file, 'w') as f_out:
  for key, value in params_out.iteritems():
    if value.find('/')==-1:
      f_out.write('%s  %s\n'%(key, value))
    else:
      f_out.write('%s "%s"\n'%(key, value))

# Get mpirun from make_kerner.macros
with open('make_kerner.macros') as f:
  for line in f:
    if line.strip()!='':
      key = line.split()[0]
      if key=='MPIRUN':
        mpirun_cmd = line.split()[2]

# Make kerner code 
subprocess.check_call('make -sj', shell=True)

# Copy code files into run_dir, tar it and delete it.
# A bit clumsy, but ensures that the internal path is Code/*.f90 etc.
code_dir = os.path.join(run_dir, 'Code')
archive_name = os.path.join(run_dir, 'Code')
os.mkdir(code_dir)
for f90_file in glob.glob('*.f90'):
  shutil.copy(f90_file, code_dir)
shutil.copy('Makefile', code_dir)
shutil.copy('make_kerner.macros', code_dir)
shutil.make_archive(archive_name, 'gztar', code_dir)
shutil.rmtree(code_dir)

shutil.copy('./kerner', run_dir)

if args.queue == 'local':
  # Change dir and submit
  os.chdir(run_dir)

  #run_cmd = mpirun_cmd, ' ../kerner', ' -n %d'%args.nslaves, ' inparam']
  run_cmd = '%s -n %d ./kerner inparam 2>&1 > OUTPUT_0000 &'%(mpirun_cmd, args.nslaves + 1)
  print run_cmd
  subprocess.call(run_cmd, shell=True)

elif args.queue == 'SuperMUC':
  # Create a LoadLeveler job script for SuperMUC
  job_script = os.path.join(run_dir, 'job.cmd')
  with open(job_script, 'w') as f:
    text_out = """
# Job file automatically created by submit.py on %s
#@ output = job_$(jobid).out 
#@ error = job_$(jobid).err 
#@ job_type = MPICH
#@ network.MPI = sn_all,not_shared,us 
#@ notification=always
#@ notify_user = staehler@geophysik.uni-muenchen.de
#@ energy_policy_tag = Kerner_intel_mpi
#@ minimize_time_to_solution = yes    
#@ class = fat"""%str(datetime.datetime.now())
    f.write(text_out)
    text_out = "#@ total_tasks=%d\n"%int(args.nslaves + 1)
    f.write(text_out)
    text_out = "#@ node = %d\n"%int(args.nslaves/40)
    f.write(text_out)
    text_out = "#@ wall_clock_limit = %d:00:00\n"%args.wall_time
    f.write(text_out)
    text_out = "#@ job_name = %s\n"%args.job_name
    f.write(text_out)
    text_out = "#@ initialdir = %s\n"%os.path.realpath(run_dir)
    f.write(text_out)
    text_out = """
#@ queue
. /etc/profile
. /etc/profile.d/modules.sh
module unload mpi.ibm
module load mpi.intel
module load fortran/intel
module load netcdf
module load fftw
module load mkl
mpiexec -n %d ./kerner inparam_basic > OUTPUT_0000\n"""%int(args.nslaves + 1)
    f.write(text_out)
  print 'Submitting to SuperMUC loadleveler queue'
  #subprocess.call(['llsubmit', job_script])

elif args.queue == 'monch':
    with open(os.path.join(run_dir, 'sbatch.sh'), 'w') as f:
        text_out ="#!/bin/bash -l\n"
        text_out += "#SBATCH --ntasks=%d\n" % (args.nslaves + 1) 
        text_out += "#SBATCH --ntasks-per-node=%d\n" % int(args.nslaves/20.)
        text_out += "#SBATCH --nodes=%d\n" % int((args.nslaves+1)/int(args.nslaves/20.))
        text_out += "#SBATCH --mem-per-cpu=%d\n" % int(args.available_memory)
        text_out += "#SBATCH --time=%d:00:00\n" % args.wall_time
        text_out += "#SBATCH --partition=_compute\n"
        text_out += "#SBATCH --job-name=%s\n" % args.job_name
        text_out += "#SBATCH --output=kerner_out.o\n"
        text_out += "#SBATCH --error=kerner_err.o\n"
        text_out += "echo The current job ID is $SLURM_JOB_ID\n"
        text_out += "echo Running on $SLURM_JOB_NUM_NODES nodes\n"
        text_out += "echo Using $SLURM_NTASKS_PER_NODE tasks per node\n"
        text_out += "echo A total of $SLURM_NTASKS tasks is used\n"
        text_out += "mpirun -n %d ./kerner inparam_basic &> OUTPUT_0000\n" % (args.nslaves + 1)
        f.write(text_out)
    os.chdir(run_dir)
    run_cmd = 'sbatch sbatch.sh'
    print run_cmd
    subprocess.call(run_cmd, shell=True)


