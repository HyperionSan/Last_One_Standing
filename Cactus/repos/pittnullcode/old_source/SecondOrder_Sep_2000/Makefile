#! /usr/sbin/smake
#
FC	= cf77
F90	= f90
CC	= cc
LD	= $(F90)
MOD	= mod
SWP	= 
RM	= /bin/rm -f
MP	=
ABI	=
ISA	=
ARCH	= $(MP) $(ABI) $(ISA)
OLEVEL	= -dp -ef -O vector3 -O scalar3
FOPTS	=
COPTS	=
F90FLAGS= $(ARCH) $(OLEVEL) $(FOPTS)
FFLAGS	= $(ARCH) $(OLEVEL) $(FOPTS)
CFLAGS	= $(ARCH) $(OLEVEL) $(COPTS)
LIBPATH = #-L./lib
LIBS	= #-lgft
LDFLAGS	= $(ARCH) $(OLEVEL) -lperf
LIBFLAGS = $(LIBPATH) $(LIBS)
FIXED	= -f fixed
FREE	= -f free

PROG =	full

SRCS =	affine.f90 ascii_io.f90 boundary.f90 calcconfactnorth.f90 \
	calcconfactsouth.f90 calcnewsnorth.f90 calcnewsouth.f90 eth.f90 \
	geomview.f90 gridtranslator.f90 hdriver.f90 hio.f90 horizon.f90 \
	horizon_eth.f90 mod_coord.f90 model.f90 modnews.f90 null.f90 \
	null_analytic.f90 null_code.f90 null_coortranvars.f90 null_eth.f90 \
	null_evol.f90 null_grid.f90 null_hyper_beta.f90 null_hyper_u.f90 \
	null_hyper_w.f90 null_interp.f90 null_io.f90 null_mask.f90 \
	null_newsvarnorth.f90 null_newsvarsouth.f90 null_params.f90 \
	null_vars.f90 numservicef90.f90 odedriver.f90 rksuite_90.f90 \
	set_precision.f90 spheroid_initial.f90 dpsort.f

OBJS =	affine.o ascii_io.o boundary.o calcconfactnorth.o calcconfactsouth.o \
	calcnewsnorth.o calcnewsouth.o eth.o geomview.o gridtranslator.o \
	hdriver.o hio.o horizon.o horizon_eth.o mod_coord.o model.o modnews.o \
	null.o null_analytic.o null_code.o null_coortranvars.o null_eth.o \
	null_evol.o null_grid.o null_hyper_beta.o null_hyper_u.o \
	null_hyper_w.o null_interp.o null_io.o null_mask.o \
	null_newsvarnorth.o null_newsvarsouth.o null_params.o null_vars.o \
	numservicef90.o odedriver.o rksuite_90.o set_precision.o \
	spheroid_initial.o dpsort.o

CNTRL = *.in makemake Make* Scripts 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBFLAGS)

clean:
	rm -f $(PROG) $(OBJS) *.$(MOD)

tar:
	tar cf `basename $(PWD)`.tar $(SRCS) *.in Make*

.SUFFIXES: $(SUFFIXES) .f90 .f

.f90.o:
	$(F90) $(FREE) $(F90FLAGS) -c $<

.f.o:
	$(F90) $(FIXED) $(F90FLAGS) -c $<

affine.o: model.o set_precision.o
	$(F90) $(F90FLAGS) -c -p model.o -p set_precision.o affine.f90
boundary.o: hdriver.o horizon.o null_grid.o null_params.o null_vars.o
	$(F90) $(F90FLAGS) -c -p hdriver.o -p horizon.o -p null_grid.o -p \
		null_params.o -p null_vars.o boundary.f90
calcconfactnorth.o: null_eth.o null_grid.o null_newsvarnorth.o null_vars.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_grid.o -p \
		null_newsvarnorth.o -p null_vars.o calcconfactnorth.f90
calcconfactsouth.o: null_eth.o null_grid.o null_newsvarsouth.o null_vars.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_grid.o -p \
		null_newsvarsouth.o -p null_vars.o calcconfactsouth.f90
calcnewsnorth.o: null_eth.o null_grid.o null_newsvarnorth.o null_vars.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_grid.o -p \
		null_newsvarnorth.o -p null_vars.o calcnewsnorth.f90
calcnewsouth.o: null_eth.o null_grid.o null_newsvarsouth.o null_vars.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_grid.o -p \
		null_newsvarsouth.o -p null_vars.o calcnewsouth.f90
gridtranslator.o: affine.o horizon.o model.o set_precision.o
	$(F90) $(F90FLAGS) -c -p affine.o -p horizon.o -p model.o -p \
		set_precision.o gridtranslator.f90
hdriver.o: gridtranslator.o hio.o horizon.o model.o spheroid_initial.o
	$(F90) $(F90FLAGS) -c -p gridtranslator.o -p hio.o -p horizon.o -p \
		model.o -p spheroid_initial.o hdriver.f90
hio.o: geomview.o horizon.o model.o
	$(F90) $(F90FLAGS) -c -p geomview.o -p horizon.o -p model.o hio.f90
horizon.o: horizon_eth.o
	$(F90) $(F90FLAGS) -c -p horizon_eth.o horizon.f90
horizon_eth.o: eth.o
	$(F90) $(F90FLAGS) -c -p eth.o horizon_eth.f90
mod_coord.o: null_coortranvars.o null_eth.o null_grid.o null_interp.o \
	null_newsvarnorth.o null_newsvarsouth.o null_params.o
	$(F90) $(F90FLAGS) -c -p null_coortranvars.o -p null_eth.o -p \
		null_grid.o -p null_interp.o -p null_newsvarnorth.o -p \
		null_newsvarsouth.o -p null_params.o mod_coord.f90
modnews.o: calcconfactnorth.o calcconfactsouth.o calcnewsnorth.o \
	calcnewsouth.o null_eth.o null_grid.o null_interp.o \
	null_newsvarnorth.o null_newsvarsouth.o null_vars.o
	$(F90) $(F90FLAGS) -c -p calcconfactnorth.o -p calcconfactsouth.o -p \
		calcnewsnorth.o -p calcnewsouth.o -p null_eth.o -p \
		null_grid.o -p null_interp.o -p null_newsvarnorth.o -p \
		null_newsvarsouth.o -p null_vars.o modnews.f90
null.o: ascii_io.o boundary.o mod_coord.o modnews.o null_analytic.o \
	null_code.o null_coortranvars.o null_grid.o null_io.o \
	null_newsvarnorth.o null_newsvarsouth.o null_vars.o
	$(F90) $(F90FLAGS) -c -p ascii_io.o -p boundary.o -p mod_coord.o -p \
		modnews.o -p null_analytic.o -p null_code.o -p \
		null_coortranvars.o -p null_grid.o -p null_io.o -p \
		null_newsvarnorth.o -p null_newsvarsouth.o -p null_vars.o \
		null.f90
null_analytic.o: null_grid.o null_params.o null_vars.o
	$(F90) $(F90FLAGS) -c -p null_grid.o -p null_params.o -p null_vars.o \
		null_analytic.f90
null_code.o: null_eth.o null_evol.o null_grid.o null_hyper_beta.o \
	null_hyper_u.o null_hyper_w.o null_interp.o null_params.o null_vars.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_evol.o -p null_grid.o -p \
		null_hyper_beta.o -p null_hyper_u.o -p null_hyper_w.o -p \
		null_interp.o -p null_params.o -p null_vars.o null_code.f90
null_coortranvars.o: null_grid.o
	$(F90) $(F90FLAGS) -c -p null_grid.o null_coortranvars.f90
null_eth.o: null_grid.o
	$(F90) $(F90FLAGS) -c -p null_grid.o null_eth.f90
null_evol.o: null_eth.o null_grid.o null_mask.o null_params.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_grid.o -p null_mask.o -p \
		null_params.o null_evol.f90
null_grid.o: null_params.o
	$(F90) $(F90FLAGS) -c -p null_params.o null_grid.f90
null_hyper_beta.o: null_grid.o null_mask.o
	$(F90) $(F90FLAGS) -c -p null_grid.o -p null_mask.o \
		null_hyper_beta.f90
null_hyper_u.o: null_eth.o null_grid.o null_mask.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_grid.o -p null_mask.o \
		null_hyper_u.f90
null_hyper_w.o: null_eth.o null_grid.o null_mask.o
	$(F90) $(F90FLAGS) -c -p null_eth.o -p null_grid.o -p null_mask.o \
		null_hyper_w.f90
null_interp.o: null_grid.o
	$(F90) $(F90FLAGS) -c -p null_grid.o null_interp.f90
null_io.o: ascii_io.o null_grid.o null_params.o null_vars.o
	$(F90) $(F90FLAGS) -c -p ascii_io.o -p null_grid.o -p null_params.o \
		-p null_vars.o null_io.f90
null_newsvarnorth.o: null_grid.o
	$(F90) $(F90FLAGS) -c -p null_grid.o null_newsvarnorth.f90
null_newsvarsouth.o: null_grid.o
	$(F90) $(F90FLAGS) -c -p null_grid.o null_newsvarsouth.f90
null_vars.o: null_grid.o
	$(F90) $(F90FLAGS) -c -p null_grid.o null_vars.f90
numservicef90.o: numservicef90.o set_precision.o
	$(F90) $(F90FLAGS) -c -p numservicef90.o -p set_precision.o \
		numservicef90.f90
odedriver.o: affine.o gridtranslator.o model.o numservicef90.o rksuite_90.o \
	set_precision.o
	$(F90) $(F90FLAGS) -c -p affine.o -p gridtranslator.o -p model.o -p \
		numservicef90.o -p rksuite_90.o -p set_precision.o \
		odedriver.f90
rksuite_90.o: rksuite_90.o
	$(F90) $(F90FLAGS) -c -p rksuite_90.o rksuite_90.f90
spheroid_initial.o: gridtranslator.o horizon_eth.o model.o
	$(F90) $(F90FLAGS) -c -p gridtranslator.o -p horizon_eth.o -p model.o \
		spheroid_initial.f90
