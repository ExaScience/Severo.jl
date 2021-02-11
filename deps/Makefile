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

ifndef LAPACK_CFLAGS
LAPACK_CFLAGS=-DHAVE_F77_UNDERSCORE
endif

ifndef LAPACK_LDFLAGS
#LAPACK_FLAGS=-lmkl_rt
LAPACK_LDFLAGS=-llapack -lblas
endif

LDFLAGS=-fopenmp $(LAPACK_LDFLAGS)
CFLAGS=-O3 -march=native -fPIC $(LAPACK_CFLAGS)
CXXFLAGS=-O3 -march=native -fPIC -fopenmp

OBJECTS=neighbours.o clustering.o ModularityOptimizer.o irlba.o

.cpp.o:
	$(CXX) -std=c++17 -c $(CXXFLAGS) $(INCPATH) -o $@ $<

.c.o:
	$(CC) -Wall -ansi -std=gnu99 -c $(CFLAGS) $(INCPATH) -o $@ $<

all: $(OBJECTS)
	$(LD) -shared $(OBJECTS) $(LDFLAGS) -o libcell.$(SLIB)

clean:
	rm $(OBJECTS) libcell.$(SLIB)
