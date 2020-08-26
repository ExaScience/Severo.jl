all:
	cc -DHAVE_F77_UNDERSCORE -Wall -ansi -pedantic -std=gnu99 -g -shared irlba.c -oirlba.so -llapack -lblas
