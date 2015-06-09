libperturbator.a: core.o generated/z2c.o
	ar rcsv libperturbator.a core.o generated/*.o

core.o: core.c generated/z2c.h
	gcc -std=c99 -pedantic -Wall -Wextra -fopenmp -O3 -march=native -c core.c

generated/z2c.o: generated/z2c.h
	make -C generated

generated/z2c.h: perturbator
	-rm -r generated
	mkdir generated
	cd generated && ../perturbator

perturbator: perturbator.hs
	ghc -O2 perturbator.hs
