############################################################################
# 'Makefile for building the Potential of Mean Force code'
# Author:  Otilio Enrique Rodríguez López (2023)
# Email: oenriquerdzlpz@gmail.com
# Department: Instituto de Física, Facultad de Ciencias, UASLP.
############################################################################

CC = icc
GCC = gcc-9 #gcc working version
file_names = potential.cpp read_parameters.cpp functions.cpp
ARGS = -diag-disable=10441

all: pfm_coulomb.x

pfm_coulomb.x: $(file_names)
	$(CC) -gcc-name=$(GCC)-O2 $(file_names) -o pfm_coulomb.x $(ARGS)
	
debug: $(file_names)
	$(CC) -gcc-name=$(GCC)-O2 $(file_names) -o pfm_coulomb.x $(ARGS) -g

clean:
	rm -f pfm_coulomb.x movie.xml *.out electrolyte_movie.xyz

### Notes ###
# For icpc use -gxx-name=g++-9 version flag
