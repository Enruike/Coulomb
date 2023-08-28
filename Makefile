############################################################################
# 'Makefile for building the Potential of Mean Force code'
# Author:  Otilio Enrique Rodríguez López (2023)
# Email: oenriquerdzlpz@gmail.com
# Department: Instituto de Física, Facultad de Ciencias, UASLP.
############################################################################

CC = icc
file_names = potential.cpp read_parameters.cpp
ARGS = -diag-disable=10441

all: pfm_coulomb.x

pfm_coulomb.x: $(file_names)
	$(CC) -O2 $(file_names) -o pfm_coulomb.x $(ARGS)

clean:
	rm -f pfm_coulomb.x