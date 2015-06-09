todo
====

phase 1
* generalize number types in core.c to core_native.c
* choose appropriate number type for pixel spacing
* choose precision of reference by pixel spacing
* approx cache uses mixed number types
* switch to new reference if all pixels are glitched

phase 2
* don't use shape test
* generate nucleus code
* generate box period code
* remove dependency on mandelbrot-numerics

phase 3
* generate per-pixel perturbation code
* generate more formulas

phase 4
* implement f??i16,32,64
* implement plain, benchmark vs perturbed
* gtk3 example program

phase 5
* GPU acceleration
