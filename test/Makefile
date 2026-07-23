# Copyright (c) 2012, Johannes Feist
# licensed under the MIT open source license, see LICENSE file

TEST_EXE := test_field_values
TEST_SRC := test_field_values.f90
GNUPDEPS := $(shell grep '\.dat' tests.gnup | sed "s/\#.*//;s/^[^']*'//;s/'.*//"|sort |uniq)
INFIS    := $(wildcard *.in)
GENFIS   := $(INFIS:.in=_fields.dat) $(INFIS:.in=_fourier.dat)
PLOTLAS  := ../bin/plotlaserfield
PLOTFOU  := ../bin/plotlaserfourier

DEFAULT: test

.PHONY: DEFAULT test plot clean

test: $(TEST_EXE) laserfields_tests.pdf
	./$(TEST_EXE)

$(TEST_EXE): $(TEST_SRC)
	$(FC) $(FCFLAGS) -I../lib -I../src $< -L../lib -llaserfields -o $@

plot: laserfields_tests.pdf

laserfields_tests.pdf : tests.gnup ${GNUPDEPS}
	gnuplot $<

%_fields.dat : %.in ${PLOTLAS}
	${PLOTLAS} $< > $@

%_fourier.dat : %.in ${PLOTFOU}
	${PLOTFOU} $< > $@

clean:
	$(RM) $(GENFIS) laserfields_tests.pdf $(TEST_EXE)
