
LIB=../lib/libfftpack.a

UNAMES := $(shell uname -s)

ifeq ($(UNAMES),Linux)

#  PGI := $(shell pgf90 2>&1)

#  ifeq ($(PGI),pgf90-Warning-No files to process)

#    F90 := pgf90 -module ../lib -I../lib
#   F90 := pgf90 -r8 -module ../lib -I../lib
#    CPP := pgf90 -E

#  else

    F90 := gfortran -fdefault-real-8 -ffast-math -ffree-line-length-0 -O3 -I../lib 
    CPP := gfortran -x f95-cpp-input 

#    F90 := g95 -DG95 -g -fmod=../lib -I../lib 
#    CPP := g95 -E -DG95

#  endif

  MAKE := make
  AR := /usr/bin/ar

endif

ifeq ($(UNAMES),AIX)

  F90 := xlf -qmoddir=../lib -I../lib
# F90 := xlf -qautodbl=dbl -qmoddir=../lib -I../lib -lxlfpmt4 -lxlfpmt8
  CPP := xlf -d -qnoobject
  MAKE := gmake
  AR := /usr/bin/ar

endif

ifeq ($(UNAMES),SunOS)

    AR := /usr/ccs/bin/ar
    F90 := /opt/SUNWspro/bin/f90 -moddir=../lib -I../lib
    CPP := /opt/SUNWspro/bin/f90 -F
    MAKE := /fs/local/bin/make

endif

ifeq ($(UNAMES),IRIX64)

    AR := /usr/bin/ar
    F90 := f90 -I../lib
    CPP := f90 -E
    MAKE := /usr/local/bin/gmake

endif

ifeq ($(UNAMES),Darwin)

    AR := /usr/bin/ar
    F90 := gfortran
#    F90 := gfortran -fdefault-real-8
    CPP := gfortran -cpp
    MAKE := /usr/bin/make

endif
