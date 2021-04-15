SHELL=/bin/sh
#
# This makefile was produced by /usr/bin/fmgen at 05:40:22 PM on 06/17/93
# If it is invoked by the command line
#	make -f makefile
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable named ecosys.x.
#
# If it is invoked by the command line
#	make -f makefile ecosys.x.prof
# it will compile the fortran modules indicated by SRCS into the object
# modules indicated by OBJS and produce an executable which profiles
# named ecosys.x.prof.
#
# To remove all the objects but leave the executables use the command line
#	make -f makefile clean
#
# To remove everything but the source files use the command line
#	make -f makefile clobber
#
# To remove the source files created by /usr/bin/fmgen and this makefile
# use the command line
#	make -f makefile void
#
# The parameters SRCS and OBJS should not need to be changed.  If, however,
# you need to add a new module add the name of the source module to the
# SRCS parameter and add the name of the resulting object file to the OBJS
# parameter.  The new modules are not limited to fortran, but may be C, YACC,
# LEX, or CAL.  An explicit rule will need to be added for PASCAL modules.
#
SRCSF := $(wildcard ./*.f)
SRCSC := $(wildcard ./*.c)
OBJSF := $(SRCSF:./%.f=../bld/%.f.o)
OBJSC := $(SRCSC:./%.c=../bld/%.c.o)
# Tunable parameters
#
# CF		Name of the fortran compiling system to use
# LDFLAGS	Flags to the loader
# LIBS		List of libraries
# CMD		Name of the executable
# PROFLIB	Library needed for profiling
#
CF =		ifort
CC =            icc
#CF =gf6
#CC =gcc6
#LIBS =		-lpgsse2
# CF =		ifort
# LIBS =
CMD =		ecosys.x
#PROFLIB =	-lprof

# To perform the default compilation, use the first line
# To compile with flowtracing turned on, use the second line
# To compile giving profile additonal information, use the third line
# WARNING:  SIMULTANEOUSLY PROFILING AND FLOWTRACING IS NOT RECOMMENDED
#FFLAGS = -axP -xW -ipo -O3 -r8 -i4 -align dcommons -cpp -save
#intel setup
FFLAGS = -O2 -mp1 -r8 -i4 -align dcommons -cpp -auto-scalar -fimf-arch-consistency=true
FFLAGS = -O2 -mp1 -r8 -i4 -align dcommons -auto-scalar -fimf-arch-consistency=true
#gfortran
#FFLAGS = -O2

#FFLAGS = -O2 -check uninit -check bounds -mp1 -r8 -i4 -align dcommons -cpp -auto-scalar
#FFLAGS = -O2 -check uninit -mp1 -r8 -i4 -align dcommons -cpp -auto-scalar
#FFLAGS = -O2 -check bounds -mp1 -r8 -i4 -align dcommons -cpp -auto-scalar
#FFLAGS = -debug -fpe0 -traceback -r8 -i8 -align dcommons -cpp -save -check bounds
OPTSC  =  -O2
#OPTSC  =  -O2  -cpp -save
#FFLAGS =	 -g -r8 -align dcommons -cpp
#FFLAGS =	 -g -r8 -align dcommons -check bounds -cpp
#FFLAGS =	 -n32 -r8 -woff 2345 -r10000 -mips4 -O2 -OPT:Olimit=6000
#FFLAGS =	 -n32 -r8 -woff 2345 -r10000 -mips4 -O0 -g -DEBUG:trap_uninitialized=ON
#FFLAGS =	 f90 *.f -O2 -o ecosys.x -lmpi
LDFLAGS = 	$(FFLAGS)

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build ecosys.x.
#
all:		$(CMD)

$(OBJSF): ../bld/%.f.o : ./%.f
	$(CF)  $(FFLAGS) -c $< -o $@
$(OBJSC): ../bld/%.c.o : ./%.c
	$(CC) -c $(OPTSC) $< -o $@
#.f.o:
#	$(CF) $(FFLAGS) -c $< -o ../bld/$@
#.c.o:
#	$(CC) -c $(OPTSC) $< -o ../bld/$@

$(CMD):$(OBJSF) $(OBJSC)
	$(CF) $(LDFLAGS) -o ../bld/$(CMD) $(OBJSF) $(OBJSC) $(LIBS)
#

# Make the profiled version of the command and call it ecosys.x.prof
#
$(CMD).prof:	$(OBJSF) $(OBJSC)
	$(CF) $(LDFLAGS) -o ../bld/$(@) $(OBJSF) $(OBJSC) $(PROFLIB) $(LIBS)

clean:
	-rm -f ../bld/$(CMD) ../bld/*.o *.il

clobber:	clean
	-rm -f $(CMD) ../bld/$(CMD).prof

void:	clobber
	-rm -f $(SRCS) makefile
