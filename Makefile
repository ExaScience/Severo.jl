ANN_FLAGS=-Iann_1.1.2/include -Lann_1.1.2/lib

all:
	cc -DHAVE_F77_UNDERSCORE -Wall -ansi -pedantic -std=gnu99 -O3 -shared irlba.c -oirlba.so -llapack -lblas
	c++ -Wall -ansi -pedantic -std=c++17 -O3 $(ANN_FLAGS) -shared neighbours.cpp -oneighbours.so -lANN
