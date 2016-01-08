# perturbator -- efficient deep zooming for Mandelbrot sets
# Copyright (C) 2015,2016 Claude Heiland-Allen
# License GPL3+ http://www.gnu.org/licenses/gpl.html

libperturbator.a: perturbator.o generated/z2c.o
	ar rcsv libperturbator.a perturbator.o generated/*.o

perturbator.o: perturbator.c generated/z2c.h perturbator.h edouble.cc
	g++ -std=c++11 -pedantic -Wall -Wextra -fopenmp -O3 -march=native -c perturbator.c -I../../mandelbrot-numerics/c/include -I../../mandelbrot-numerics/c/lib -ggdb

generated/z2c.o: generated/z2c.h edouble.cc
	make -C generated

generated/z2c.h: generator
	-rm -r generated
	mkdir generated
	cd generated && ../generator

generator: generator.hs
	ghc -O2 generator.hs
