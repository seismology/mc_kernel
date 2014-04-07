PROG = rdbm

SRCS =	commpi.F90 buffer.f90 clocks.f90 fft_type.f90 filtering.f90 \
	global_parameters.f90 kdtree2.f90 kernel.f90 nc_routines.f90 rdbm.f90 \
	readfields.f90 resampling_type.f90 simple_routines.f90 \
	sourceandreceiver.f90 test_resampling.f90

OBJS =	commpi.o buffer.o clocks.o fft_type.o filtering.o global_parameters.o \
	kdtree2.o kernel.o nc_routines.o rdbm.o readfields.o \
	resampling_type.o simple_routines.o sourceandreceiver.o \

#Example to include specific netcdf libraries: 
LIBS = -lm -llapack -lfftw3 -lfftw3f -lnetcdff 

# set unc to compile with netcdf: 
#F90FLAGS = -Dunc 
CC = gcc
CFLAGS = -O3 -DF_UNDERSCORE
F90 = mpif90 -O3  -fbacktrace -g 
FC = gfortran -O3 -fbacktrace -g 
 
# to include local built of netcdf you might want to use sth like this: 
INCLUDE = -I /usr/include 


# cancel m2c implicit rule 
%.o : %.mod 
 

all: $(PROG) test_resampling

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

test_resampling: $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.M *.mod *.d *.il core 

.SUFFIXES: $(SUFFIXES) .f90 .F90

.f90.o:
	$(F90) $(F90FLAGS) -c $(INCLUDE) $<

.F90.o:
	$(F90) $(F90FLAGS) -c $(INCLUDE) $<

buffer.o: commpi.o global_parameters.o Makefile 
clocks.o: global_parameters.o Makefile 
fft_type.o: commpi.o global_parameters.o simple_routines.o  Makefile 
filtering.o: commpi.o fft_type.o global_parameters.o simple_routines.o Makefile 
global_parameters.o:  Makefile 
kdtree2.o: kdtree2.o Makefile 
kernel.o: commpi.o fft_type.o filtering.o global_parameters.o Makefile 
nc_routines.o: commpi.o global_parameters.o Makefile 
rdbm.o: clocks.o commpi.o fft_type.o global_parameters.o readfields.o \
	resampling_type.o sourceandreceiver.o Makefile 
readfields.o: buffer.o clocks.o commpi.o global_parameters.o kdtree2.o \
	nc_routines.o simple_routines.o sourceandreceiver.o Makefile 
resampling_type.o: commpi.o fft_type.o global_parameters.o Makefile 
simple_routines.o: global_parameters.o Makefile 
sourceandreceiver.o: global_parameters.o kernel.o Makefile 
test_resampling.o: commpi.o global_parameters.o resampling_type.o Makefile 
commpi.F90: global_parameters.o Makefile 
