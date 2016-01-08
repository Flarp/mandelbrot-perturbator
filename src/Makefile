# perturbator -- efficient deep zooming for Mandelbrot sets
# Copyright (C) 2015,2016 Claude Heiland-Allen
# License GPL3+ http://www.gnu.org/licenses/gpl.html

libperturbator.a: core.o generated/z2c.o
	ar rcsv libperturbator.a core.o generated/*.o

core.o: core.c core_native.c generated/z2c.h perturbator.h edouble.cc
	g++ -std=c++11 -pedantic -Wall -Wextra -fopenmp -O3 -march=native -c core.c -I/home/claude/code/borok/mandelbrot-numerics/c/include -I/home/claude/code/borok/mandelbrot-numerics/c/lib -ggdb

generated/z2c.o: generated/z2c.h edouble.cc
	make -C generated

generated/z2c.h: perturbator
	-rm -r generated
	mkdir generated
	cd generated && ../perturbator

perturbator: perturbator.hs
	ghc -O2 perturbator.hs
