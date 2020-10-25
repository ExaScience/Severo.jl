ifeq ($(OS), Windows_NT)
	CC=gcc
	CXX=g++
	LD=g++
	SLIB=dll
else
	UNAME := $(shell uname)
	ifeq ($(UNAME), Darwin)
		CC=gcc
		CXX=g++
		LD=g++
		SLIB=dylib
	else
		CC=gcc
		CXX=g++
		LD=g++
		SLIB=so
	endif
endif

#LAPACK_FLAGS=-fopenmp -lmkl_rt
#LAPACK_FLAGS=-L/opt/software/OpenBLAS/0.3.7-GCC-8.3.0/lib -lopenblas -Wl,-rpath=/opt/software/OpenBLAS/0.3.7-GCC-8.3.0/lib
LAPACK_FLAGS=-llapack -lblas

CFLAGS=-DHAVE_F77_UNDERSCORE -O3 -march=native -fPIC
CXXFLAGS=-O3 -march=native -fPIC

OBJECTS=neighbours.o clustering.o ModularityOptimizer.o irlba.o

.cpp.o:
	$(CXX) -std=c++17 -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -Wall -ansi -std=gnu99 -c $(CFLAGS) $(INCPATH) -o $@ $<

all: $(OBJECTS)
	$(LD) -shared $(OBJECTS) -o libcell.$(SLIB)

clean:
	rm $(OBJECTS) libcell.$(SLIB)
