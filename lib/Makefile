# perturbator -- efficient deep zooming for Mandelbrot sets
# Copyright (C) 2015,2016 Claude Heiland-Allen
# License GPL3+ http://www.gnu.org/licenses/gpl.html

libperturbator.a: perturbator.o
	ar rcsv libperturbator.a perturbator.o

perturbator.o: perturbator.c perturbator.h edouble.cc z2c.c
	g++ -std=c++11 -pedantic -Wall -Wextra -fopenmp -fPIC -O3 -march=native -c perturbator.c -I../../mandelbrot-numerics/c/include -I../../mandelbrot-numerics/c/lib -ggdb

.PHONY: clean
clean:
	-rm -r libperturbator.a perturbator.o
