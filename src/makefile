# !! Define Environment variable FARGO_ARCH to reflect your
#architecture for ex.: setenv FARGO_ARCH LINUX (or export
#FARGO_ARCH=LINUX depending on your shell) 

# Clement's Desktop computer at IRAP
# Opteron platform FARGO_ARCH must be set to CIRAP
CC_CIRAP = gcc
OPT_CIRAP = -O3 -m64 -ffast-math
OPTSEQ_CIRAP =
PARAOPT_CIRAP =
PARACC_CIRAP = mpicc

# Titan cluster at IRAP
# Opteron platform FARGO_ARCH must be set to TITAN
CC_TITAN = icc
OPT_TITAN = -O3 -m64 -Wno-unknown-pragmas -Wall
OPTSEQ_TITAN =
PARAOPT_TITAN = 
PARACC_TITAN = mpicc

# CALMIP regional cluster at University of Toulouse
# Opteron platform FARGO_ARCH must be set to CALMIP
CC_CALMIP = icc
OPT_CALMIP = -O3 -m64 -Wno-unknown-pragmas -Wall
OPTSEQ_CALMIP =
PARAOPT_CALMIP = 
PARACC_CALMIP = mpiicc

# Macintosh MAC OS X platform (with gcc) FARGO_ARCH must be set to
#MACOSX
CC_MACOSX  = cc
OPT_MACOSX = -O3 -m64 -ffast-math 
OPTSEQ_MACOSX = 
PARAOPT_MACOSX =
PARACC_MACOSX = mpicc

# Generic options valid for any platfom which has 'cc' in the path
#These are used if you leave FARGO_ARCH undefined
CC_  = cc
OPT_ = -O3
OPTSEQ_ = 
PARAOPT_ =
PARACC_ = mpicc

#
#
#
#--------------------No Changes Needed after this line------------------------
#
#
#
SHELL		=  /bin/sh

MAINOBJ         = LowTasks.o SideEuler.o Output.o Init.o main.o Theo.o\
		  Interpret.o SourceEuler.o TransportEuler.o Stockholm.o\
		  Planet.o RungeKunta.o Viscosity.o Psys.o Force.o var.o\
		  Pframeforce.o split.o merge.o commbound.o fpe.o rebin.o\
		  sgmain.o sginit.o sgdens.o sgkernel.o sgacc.o sgzero.o\
		  sgupdate.o sgsysinit.o axilib.o aniso.o turb.o\
		  Dsys.o DustUpdate.o \

MPIDUMMY	= mpi_dummy.o
FFTWDUMMY	= fftw_dummy.o

COMP        = $(CC_$(FARGO_ARCH))
PARACOMP    = $(PARACC_$(FARGO_ARCH))
OPT         = $(OPT_$(FARGO_ARCH))
OPTSEQ      = $(OPTSEQ_$(FARGO_ARCH))
PARAOPT     = $(PARAOPT_$(FARGO_ARCH)) -D_PARALLEL -I$(MPI_PREFIX)/include
FFTWOPT	    = -D_FFTW -I$(FFTW_PREFIX)/include
LIBS        = -lm
PARALIBS    = -L$(MPI_PREFIX)/lib -lmpi
FFTWLIBS    = -L$(FFTW_PREFIX)/lib -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw
AUTOINCL    = param.h param_noex.h global_ex.h sg_ex.h

include	.config

EXENAME        = ../fargo

ARCHIVE		= $(EXENAME:../%=%.tar)
ARCHIVECOMP	= $(EXENAME:../%=%.tar.gz)

SRC = *.c
INCLUDE = *.h

ifeq ($(BUILD),parallelfftw)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS) $(FFTWLIBS)
OPTIONS		= $(OPT) $(PARAOPT) $(FFTWOPT)
OBJ		= $(MAINOBJ)
endif
ifeq ($(BUILD),parallel)
COMPILER	= $(PARACOMP)
LIBRARIES	= $(LIBS) $(PARALIBS)
OPTIONS		= $(OPT) $(PARAOPT)
OBJ		= $(MAINOBJ) $(FFTWDUMMY)
endif
ifeq ($(BUILD),sequential)
COMPILER	= $(COMP)
LIBRARIES	= $(LIBS)
OPTIONS		= $(OPT) $(OPTSEQ)
OBJ		= $(MAINOBJ) $(MPIDUMMY) $(FFTWDUMMY)
endif

all: conditionalrebuild $(AUTOINCL) $(OBJ) $(EXENAME) archive
	@echo "" 
	@echo ""
	@echo "      NOTE"
	@echo ""
ifeq ($(BUILD),parallelfftw)
	@echo "This built is PARALLEL (MPI) and uses FFTW librairies"
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo "gmake BUILD=parallel"
	@echo "gmake BUILD=sequential"
endif
ifeq ($(BUILD),parallel)
	@echo "This built is PARALLEL (MPI)"
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo "gmake BUILD=parallelfftw"
	@echo "gmake BUILD=sequential"
endif
ifeq ($(BUILD),sequential)
	@echo "This built is SEQUENTIAL"
	@echo "If you want to change this,"
	@echo "then you need to issue:"
	@echo "gmake BUILD=parallelfftw"
	@echo "gmake BUILD=parallel"
endif
	@echo ""

$(EXENAME): $(OBJ)
	$(COMPILER) $(OBJ) $(OPTIONS) -o $(EXENAME) $(LIBRARIES)

.PHONY: conditionalrebuild
ifneq ($(BUILD),$(OLDBUILD))
conditionalrebuild: clean
	@echo "OLDBUILD = $(BUILD)" > .config
	@echo "BUILD = $(BUILD)" >> .config
else
conditionalrebuild:
endif

.oldconfig:
.config:

archive : $(SRC) $(INCL) makefile varparser.pl	
	@echo "Creating ../source.tar.bz2"
	@tar cf ../source.tar *.c
	@tar rf ../source.tar *.h
	@tar rf ../source.tar makefile
	@tar rf ../source.tar varparser.pl
	@bzip2 -9 -f ../source.tar

para:
	@gmake BUILD=parallel
parafftw:
	@gmake BUILD=parallelfftw
seq:
	@gmake BUILD=sequential

$(AUTOINCL) : var.c global.h sg.h makefile varparser.pl
	@./varparser.pl

$(OBJ): mp.h fondam.h param.h param_noex.h types.h makefile

.PHONY: clean mrproper package

mrproper:
	rm -f *.o *~ *.s *.il $(AUTOINCL) $(EXENAME) ../core.*\
	*.tex *.dvi *.pdf *.ps *.log *.aux *.lint $(ARCHIVE)\
	$(ARCHIVECOMP)

clean:
	rm -f *.o *~ *.s *.il

package: $(ARCHIVECOMP)

# The Doxyfile provided here is relatively general and should allow
# you to build your own documentation if you modify the source. You
# should have the required executable files in your path
# (doxygen, dot, latex, epstopdf, dvips, gs). On Mac Os X I made soft
# links in /usr/bin to these executables (either in /sw/bin or
# /Applications/Doxygen.app/Content/Resources).
doc:
	doxygen Doxyfile

release:
	@echo "Creating archive fargoadsg.tar.gz for release"
	@cd ../; tar -c\
         --exclude src/adimvalue.c\
         -f fargoadsg.tar src/*.c
	@cd ..; tar rf fargoadsg.tar src/makefile
	@cd ..; tar rf fargoadsg.tar src/*.pl
	@cd ..; tar rf fargoadsg.tar src/*.h
	@cd ..; tar rf fargoadsg.tar in/*.par
	@cd ..; tar rf fargoadsg.tar in/*.cfg
	@cd ..; tar rf fargoadsg.tar idl/*.pro
	@cd ..; tar rf fargoadsg.tar out1/
	@cd ..; gzip -9 fargoadsg.tar

.c.o  :
	$(COMPILER) $*.c -c $(OPTIONS)
