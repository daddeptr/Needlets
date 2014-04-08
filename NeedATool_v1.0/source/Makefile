Healpix_LIB = /Users/dpietrob/work/cmbcodes/Healpix_2.13a/lib_gfortran

Healpix_INCLUDE = /Users/dpietrob/work/cmbcodes/Healpix_2.13a/include_gfortran

CFITSIO = /usr/local/lib

F90C    = gfortran
FFLAGS  = -O -DGFORTRAN -fno-second-underscore -m64 -fopenmp -ffree-line-length-none #$(CHECKFLAGS)
CHECKFLAGS = -fbounds-check -g -fbacktrace

LIBRARIES = -L$(Healpix_LIB) -L$(CFITSIO) -lhealpix -lcfitsio -framework veclib

LAPACKL =

SNAME = ../synneed

ANAME = ../ananeed

F90FLAGS = -O3 -I$(Healpix_INCLUDE) $(FFLAGS) 

SOBJFILES = needlets_mod10.o  synneed_dr2.o

AOBJFILES = needlets_mod10.o  ananeed_dr2.o

NeedAToolLIB = libNeedATool.a

default: synneed ananeed

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90


synneed: $(SOBJFILES) 	
	$(F90C) $(F90FLAGS) $(SOBJFILES) $(LIBRARIES) -o $(SNAME)

ananeed: $(AOBJFILES) 	
	$(F90C) $(F90FLAGS) $(AOBJFILES) $(LIBRARIES) -o $(ANAME)

natlib: needlets_mod10.o
	ar -r $(NeedAToolLIB) $?

all: synneed ananeed natlib

clean:
	rm -f *.o *.mod *.d *.pc *.obj core
tidy:
	$(clean) $(SNAME) $(ANAME) $(NeedAToolLIB)
