#LAPACK_FLAGS=-fopenmp -lmkl_rt
#LAPACK_FLAGS=-L/opt/software/OpenBLAS/0.3.7-GCC-8.3.0/lib -lopenblas -Wl,-rpath=/opt/software/OpenBLAS/0.3.7-GCC-8.3.0/lib
LAPACK_FLAGS=-llapack -lblas
CXX_FLAGS=-O3 -march=native -fPIC

all:
	cc -DHAVE_F77_UNDERSCORE -Wall -ansi -pedantic -std=gnu99 $(CXX_FLAGS) -shared irlba.c -oirlba.so $(LAPACK_FLAGS)
	c++ -Wall -pedantic -std=c++17 $(CXX_FLAGS) -shared neighbours.cpp -oneighbours.so
	c++ -std=c++17 $(CXX_FLAGS) -shared clustering.cpp ModularityOptimizer.cpp -oclustering.so
