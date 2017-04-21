# perturbator -- efficient deep zooming for Mandelbrot sets
# Copyright (C) 2015,2016,2017 Claude Heiland-Allen
# License GPL3+ http://www.gnu.org/licenses/gpl.html

libperturbator.a: perturbator.o
	ar rcsv libperturbator.a perturbator.o

perturbator.o: perturbator.cpp perturbator.h edouble.cc z2c.c
	em++ -std=c++14 -pedantic -Wall -Wextra -fopenmp -fPIC -O3 -march=native -c perturbator.cpp -I../../mandelbrot-numerics/c/include -I../../mandelbrot-numerics/c/lib -I/usr/local/include

edouble: edouble.cc
	em++ -std=c++14 -pedantic -Wall -Wextra -fopenmp -fPIC -O3 -march=native -o edouble -DEDOUBLE_STANDALONE edouble.cc -I/usr/local/include

.PHONY: clean
clean:
	-rm -r libperturbator.a perturbator.o
