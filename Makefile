ANN_FLAGS=-Iann_1.1.2/include -Lann_1.1.2/lib
LAPACK_FLAGS=-fopenmp -lmkl_rt
#LAPACK_FLAGS=-L/opt/software/OpenBLAS/0.3.7-GCC-8.3.0/lib -lopenblas -Wl,-rpath=/opt/software/OpenBLAS/0.3.7-GCC-8.3.0/lib
#LAPACK_FLAGS=-llapack -lblas

all:
	cc -DHAVE_F77_UNDERSCORE -Wall -ansi -pedantic -std=gnu99 -O3 -fPIC -shared irlba.c -oirlba.so -march=native $(LAPACK_FLAGS)
	c++ -Wall -pedantic -std=c++17 -O3 -fopenmp -fPIC $(ANN_FLAGS) -shared neighbours.cpp -oneighbours.so -lANN
	c++ -std=c++17 -O3 -fopenmp -fPIC -shared clustering.cpp ModularityOptimizer.cpp -oclustering.so
