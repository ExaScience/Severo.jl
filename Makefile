ANN_FLAGS=-Iann_1.1.2/include -Lann_1.1.2/lib

all:
	cc -DHAVE_F77_UNDERSCORE -Wall -ansi -pedantic -std=gnu99 -O3 -fPIC -shared irlba.c -oirlba.so -llapack -lblas
	c++ -Wall -pedantic -std=c++17 -O3 -fopenmp -fPIC $(ANN_FLAGS) -shared neighbours.cpp -oneighbours.so -lANN
	c++ -std=c++17 -O3 -fopenmp -fPIC -shared clustering.cpp ModularityOptimizer.cpp -oclustering.so
