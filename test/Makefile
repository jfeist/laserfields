# Copyright (c) 2012, Johannes Feist
# licensed under the MIT open source license, see LICENSE file

GNUPDEPS:=$(shell grep '\.dat' tests.gnup | sed "s/\#.*//;s/^[^']*'//;s/'.*//"|sort |uniq)
INFIS   :=${wildcard *.in}
GENFIS  :=${INFIS:.in=_fields.dat} ${INFIS:.in=_fourier.dat}
PLOTLAS :=../bin/plotlaserfield
PLOTFOU :=../bin/plotlaserfourier

laserfields_tests.pdf : tests.gnup ${GNUPDEPS}
	gnuplot $<

%_fields.dat : %.in ${PLOTLAS}
	${PLOTLAS} $< > $@

%_fourier.dat : %.in ${PLOTFOU}
	${PLOTFOU} $< > $@

clean:
	${RM} ${GENFIS} laserfields_tests.pdf
