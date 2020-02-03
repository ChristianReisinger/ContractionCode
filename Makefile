# ##########


# Makefile
CXX=g++
CXXFLAGS=-O3 -DF_
INCLUDE = -I./Lime/lime/include/
LIBS = ./Lime/lime/lib/liblime.a -llapack -lblas -lm -lgfortran

all: o_files

# ##########

o_files: contract_twopoint.o DML_crc32.o dml.o fields.o fuzz.o io.o io_utils.o propagator_io.o Q_phi.o ranlux.o ranlxd.o ranlxs.o smearing_techniques.o Wilson_loops.o

# ##########

contract_twopoint.o: contract_twopoint.cc contract_twopoint.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

DML_crc32.o: DML_crc32.cc
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

dml.o: dml.cc dml.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

fields.o: fields.cc fields.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

fuzz.o: fuzz.cc fuzz.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

io.o: io.cc io.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

io_utils.o: io_utils.cc io_utils.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

# gauge_io.o: gauge_io.cc gauge_io.hh
# 	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

propagator_io.o: propagator_io.cc propagator_io.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

Q_phi.o: Q_phi.cc Q_phi.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

ranlux.o: ranlux.cc ranlux.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o ranlux.o ranlux.cc

ranlxd.o: ranlxd.h ranlxd.c
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o ranlxd.o ranlxd.c

ranlxs.o: ranlxs.h ranlxs.c
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o ranlxs.o ranlxs.c

smearing_techniques.o: smearing_techniques.cc smearing_techniques.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $< -fopenmp

Wilson_loops.o: Wilson_loops.cc Wilson_loops.hh
	${CXX} ${INCLUDE} ${CXXFLAGS} -c -o $@ $<

# ##########

clean:
	rm -f *~ *.o
	
