COMPILER = g++
CORES = $(shell grep processor /proc/cpuinfo | wc -l)
GCC_FLAGS = -Wall -W -ansi -pedantic -Wcast-align -Wcast-qual -Wpointer-arith -Wundef -Wconversion -g -O3

SOURCES = main.cc
OBJECTS = $(foreach file,$(SOURCES),$(subst .cc,.o,$(file)))
BINARY = remesher

GCC_INCL = -I.
GCC_LINK = -Llibremesh -lremesh -lpthread -static-libstdc++
GSL_LIBS = $(shell pkg-config --libs gsl)

all: libraries ${OBJECTS}
	${COMPILER} -o ${BINARY} ${OBJECTS} ${GCC_LINK} ${GSL_LIBS}
	rm -f *.o
	rm -f libremesh/*.o
	rm -f libremesh/*.a

multicore:
	$(MAKE) -j${CORES} all

depend:
	${COMPILER} ${GCC_INCL} -MM ${SOURCES} > Makefile.dep

%.o: %.cc
	${COMPILER} -c -o $@ $< ${GCC_FLAGS} ${GCC_INCL}

clean: FORCE
	rm -f ${BINARY} ${OBJECTS}

libraries: FORCE
	$(MAKE) -C libremesh

FORCE:

include Makefile.dep

