# Copyright (c) 2012, Johannes Feist
# licensed under the MIT open source license, see LICENSE file

.SUFFIXES:
.SUFFIXES: .f90 .o .mod .a

# fortran compiler to use
FC=gfortran
# flags to compile with
FCFLAGS=-Wall -O3 -march=native
# flag that tells the compiler where to put the produced module (.mod) files
MODLOCFLAG=-J

# for intel fortran:
ifeq (${USEIFORT},yes)
  FC=ifort
  FCFLAGS=-warn all -O3 -xHOST
  MODLOCFLAG=-module
endif

LIB   :=lib/liblaserfields.a
LIBMOD:=lib/laserfields.mod
SRCS:=$(wildcard src/*.f90)
PROGSRCS:=$(wildcard programs/*.f90)
PROGOBJS:=${PROGSRCS:.f90=.o}
PROGS:=${PROGSRCS:programs/%.f90=bin/%}

DEFAULT: lib
.PHONY : lib progs clean DEFAULT
lib: ${LIB} ${LIBMOD}
progs: ${PROGS}

test: progs
	${MAKE} -C test

clean:
	${RM} ${LIB} ${LIBMOD} src/*.o src/*.mod programs/*.o bin/*
	${MAKE} -C test clean

${LIB}    : ${SRCS:.f90=.o}
${LIBMOD} : ${LIBMOD:lib/%=src/%}
	cp $< $@

${PROGOBJS} : INCFLAGS=-Ilib
${PROGOBJS} : ${LIBMOD}
%.o: %.f90
	$(FC) $(FCFLAGS) $(INCFLAGS) $(MODLOCFLAG) $(dir $@) -c $< -o $@

%.a :
	$(AR) $(ARFLAGS) $@ $^

bin/% : LDFLAGS+=-Llib -llaserfields
bin/% : programs/%.o ${LIB}
	$(FC) $(FCFLAGS) $< $(LDFLAGS) -o $@

# this tells make that it can make a .mod file by making the associated .o file
# $(NOOP) is not defined and so does nothing, but it's not an empty rule, which
# would mislead make
%.mod: %.o
	$(NOOP)

### dependencies
src/atomic_units.o : src/nrtype.mod
src/laserfields.o : src/laserfields_fileops.mod src/laserfields_module.mod
src/laserfields_fileops.o : src/laserfields_module.mod src/misc_fileops.mod src/nrtype.mod
src/laserfields_module.o : src/atomic_units.mod src/faddeeva.mod src/misc_fileops.mod src/laserfields_miscfuncs.mod src/nrtype.mod
src/misc_fileops.o : src/nrtype.mod
src/laserfields_miscfuncs.o : src/nrtype.mod
