#! /usr/sbin/smake
#
FC      = gfortran
F90     = gfortran
CC      = gcc
LD      = $(F90)
MOD     = mod
SWP     =
RM      = /bin/rm -f
MP      =
ABI     =
ISA     =
ARCH    = $(MP) $(ABI) $(ISA)
OLEVEL  = -O3
FOPTS   = -Fixed
F90OPTS = -O3
COPTS   =
F90FLAGS= $(ARCH) $(OLEVEL) $(F90OPTS)
FFLAGS  = $(ARCH) $(OLEVEL) $(FOPTS)
CFLAGS  = $(ARCH) $(OLEVEL) $(COPTS)
PROF    =

PROG =	factor

SRCS =	affine.f90 ascii_io.f90 axiio.f90 boundary.f90 cauchydump.f90 \
	checkpoint.f90 checkpoint_defs.f90 estimateTX.f90 evolve2D.f90 \
	geomview.f90 gridtranslator.f90 hdriver.f90 hio.f90 \
	horizon.f90 horizon_eth.f90 model.f90 news2.f90 null.f90 \
	null_analytic.f90 null_cfl_test.f90 null_code.f90 null_eth.f90 \
	null_evol.f90 null_grid.f90 null_hyper_beta.f90 null_hyper_u.f90 \
	null_hyper_w.f90 null_interp.f90 null_interp2.f90 null_io.f90 \
	null_mask.f90 null_params.f90 null_vars.f90 numservicef90.f90 \
	odedriver.f90 particle.f90 pointdump.f90 pureschwarz.f90 \
	rksuite_90.f90 set_precision.f90 spheroid_initial.f90 \
	StereoIntegration.f90 dpsort.f rksuite.f TS_et_d.f

OBJS =	affine.o ascii_io.o axiio.o boundary.o cauchydump.o checkpoint.o \
	checkpoint_defs.o estimateTX.o evolve2D.o geomview.o \
	gridtranslator.o hdriver.o hio.o horizon.o horizon_eth.o model.o \
	news2.o null.o null_analytic.o null_cfl_test.o null_code.o null_eth.o \
	null_evol.o null_grid.o null_hyper_beta.o null_hyper_u.o \
	null_hyper_w.o null_interp.o null_interp2.o null_io.o null_mask.o \
	null_params.o null_vars.o numservicef90.o odedriver.o particle.o \
	pointdump.o pureschwarz.o rksuite_90.o set_precision.o \
	spheroid_initial.o StereoIntegration.o dpsort.o rksuite.o TS_et_d.o

CNTRL = *.in makemake Make* Scripts 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBFLAGS)

clean:
	rm -f $(PROG) $(OBJS) *.$(MOD)

tar:
	tar cf `basename $(PWD)`.tar $(SRCS) *.in make* Make*
	gzip -f `basename $(PWD)`.tar

zip:
	zip -r `basename $(PWD)`.zip $(SRCS) *.in Make*

.SUFFIXES: $(SUFFIXES) .f90 .f .o

.f90.o:
	$(F90) $(FREE) $(F90FLAGS) -c $<

.f.o:
	$(FC) $(FIXED) $(FFLAGS) -c $<

affine.o: model.o set_precision.o
boundary.o: ascii_io.o horizon.o null_eth.o null_grid.o null_params.o \
	null_vars.o pointdump.o
cauchydump.o: checkpoint_defs.o null_grid.o null_interp.o null_params.o \
	null_vars.o
checkpoint.o: cauchydump.o checkpoint_defs.o horizon.o news2.o null_code.o \
	null_grid.o null_params.o null_vars.o
estimateTX.o: affine.o numservicef90.o
evolve2D.o: affine.o axiio.o gridtranslator.o horizon.o null_grid.o \
	null_params.o numservicef90.o set_precision.o
gridtranslator.o: affine.o model.o null_grid.o set_precision.o
hdriver.o: checkpoint.o checkpoint_defs.o evolve2D.o gridtranslator.o hio.o \
	horizon.o model.o null_params.o spheroid_initial.o
hio.o: geomview.o horizon.o model.o null_grid.o null_params.o
horizon.o: axiio.o gridtranslator.o horizon_eth.o model.o null_grid.o \
	particle.o
horizon_eth.o: null_eth.o null_interp.o
model.o: null_grid.o
news2.o: ascii_io.o checkpoint_defs.o horizon_eth.o null_grid.o null_interp.o \
	null_interp2.o null_params.o null_vars.o
null.o: StereoIntegration.o ascii_io.o boundary.o cauchydump.o checkpoint.o \
	hdriver.o model.o news2.o null_analytic.o null_cfl_test.o \
	null_code.o null_grid.o null_io.o null_params.o null_vars.o \
	particle.o pointdump.o pureschwarz.o
null_analytic.o: model.o null_grid.o null_params.o null_vars.o particle.o
null_cfl_test.o: null_eth.o null_grid.o
null_code.o: boundary.o null_cfl_test.o null_evol.o null_grid.o \
	null_hyper_beta.o null_hyper_u.o null_hyper_w.o null_interp.o \
	null_params.o null_vars.o particle.o
null_eth.o: null_grid.o
null_evol.o: null_eth.o null_grid.o null_mask.o null_params.o particle.o
null_grid.o: null_params.o
null_hyper_beta.o: null_eth.o null_grid.o null_mask.o particle.o
null_hyper_u.o: null_eth.o null_grid.o null_mask.o particle.o
null_hyper_w.o: null_eth.o null_grid.o null_mask.o particle.o
null_interp.o: null_grid.o
null_interp2.o: null_grid.o
null_io.o: ascii_io.o null_grid.o null_params.o null_vars.o
null_vars.o: null_grid.o
numservicef90.o: numservicef90.o set_precision.o
odedriver.o: affine.o gridtranslator.o model.o rksuite_90.o set_precision.o
particle.o: null_eth.o null_grid.o null_vars.o
pointdump.o: null_grid.o null_params.o
pureschwarz.o: null_grid.o null_params.o null_vars.o
rksuite_90.o: rksuite_90.o
spheroid_initial.o: axiio.o evolve2D.o gridtranslator.o horizon_eth.o model.o \
	numservicef90.o set_precision.o
