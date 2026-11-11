#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

include $(top_srcdir)/src/mk/common-libonly.am

noinst_LIBRARIES = lib.a

lib_a_SOURCES = \
	constants.fpp \
	mpi_parallel.fpp \
	messages.fpp \
	random.fpp \
	math_libs.fpp \
	string.fpp \
	table_libs.fpp \
	timers.fpp \
	ffte_fft235.fpp \
	ffte_kernel.fpp \
	ffte_zfft1d.fpp \
	fft3d_1dalltoall.fpp \
	fft3d_2dalltoall.fpp \
	fft3d_slab.fpp \
	fileio.fpp \
	fileio_control.fpp \
	fileio_data.fpp \
	fileio_data_.c \
	fileio_pdb.fpp \
	fileio_crd.fpp \
	fileio_top.fpp \
	fileio_par.fpp \
	fileio_str.fpp \
	fileio_psf.fpp \
	fileio_rst.fpp \
	fileio_prmtop.fpp \
	fileio_ambcrd.fpp \
	fileio_gropp.fpp \
	fileio_grotop.fpp \
	fileio_grocrd.fpp \
	fileio_localres.fpp \
	hardwareinfo.fpp \
	dSFMT/dSFMT.c \
	dSFMT/dSFMT_if.c \
	dSFMT/dSFMT-common.h \
	dSFMT/dSFMT-params.h \
	dSFMT/dSFMT-params11213.h \
	dSFMT/dSFMT-params1279.h \
	dSFMT/dSFMT-params132049.h \
	dSFMT/dSFMT-params19937.h \
	dSFMT/dSFMT-params216091.h \
	dSFMT/dSFMT-params2203.h \
	dSFMT/dSFMT-params4253.h \
	dSFMT/dSFMT-params44497.h \
	dSFMT/dSFMT-params521.h \
	dSFMT/dSFMT-params86243.h \
	dSFMT/dSFMT.h \
	molecules_str.fpp \
	select_atoms_str.fpp \
	fitting_str.fpp \
	molecules.fpp \
	dihedral_libs.fpp \
	select_lexer.fpp \
	select_parser.fpp \
	select_contacts.fpp \
	select_atoms.fpp \
	select.fpp \
	fitting.fpp

EXTRA_DIST = \
	Makefile.depends 

SUFFIXES = .fpp .c
DEPENDS	= Makefile.depends

CLEANFILES = $(noinst_LIBRARIES) $(lib_a_OBJECTS) *.f90 *~ *.mod

.fpp.o:
if USE_NVHPC
	$(FPP2F90_NVHPC) $*.fpp > $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
else
if USEKCOMP
	cp $*.fpp $*.f90
	$(FPP2F90) $*.f90
	mv $*.cpp.f90 $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
else
	$(FPP2F90) $*.fpp $*.f90
	$(F90COMPILE) $(FCFLAGS) -c $*.f90
endif
endif

ffte_kernel.o: FCFLAGS=$(FCFLAGS_FFTE_KERNEL_FPP__PGI__SP)

lib.a : $(lib_a_OBJECTS)
	$(AR) $(ARFLAGS) cr $(noinst_LIBRARIES) $(lib_a_OBJECTS)

depend: clean_depend
	python3 ../../fortdep.py *.fpp > $(DEPENDS)

clean_depend:
	rm -f $(DEPENDS)

-include $(DEPENDS)
