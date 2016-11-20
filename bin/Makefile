# perturbator -- efficient deep zooming for Mandelbrot sets
# Copyright (C) 2015,2016 Claude Heiland-Allen
# License GPL3+ http://www.gnu.org/licenses/gpl.html

perturbator-glfw3: glfw3.c ../lib/libperturbator.a
	gcc -std=c99 -pedantic -Wall -Wextra -fopenmp -pthread -O3 -march=native -o perturbator-glfw3 glfw3.c -I../lib -L../lib -I$(HOME)/opt/include -L$(HOME)/opt/lib -lperturbator -lmandelbrot-numerics -lGL -lGLEW -lglfw -lmpc -lmpfr -lgmp -lm -lstdc++ -ggdb

perturbator-offline: offline.c ../lib/libperturbator.a
	gcc -std=c99 -pedantic -Wall -Wextra -fopenmp -pthread -O3 -march=native -o perturbator-offline offline.c -I../lib -L../lib -I$(HOME)/opt/include -L$(HOME)/opt/lib -lperturbator -lmandelbrot-numerics -lmpc -lmpfr -lgmp -lm -lstdc++ -ggdb

,PHONY: clean
clean:
	-rm perturbator-glfw3 perturbator-offline
