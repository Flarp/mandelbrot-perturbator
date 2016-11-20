// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015,2016 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#define _POSIX_C_SOURCE 199309L

#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include <mpfr.h>
#include <mpc.h>

#include "perturbator.h"

static inline int max(int a, int b) {
  return a > b ? a : b;
}

static int envi(const char *name, int def) {
  const char *e = getenv(name);
  if (e) {
    return atoi(e);
  } else {
    return def;
  }
}

static double envd(const char *name, double def) {
  const char *e = getenv(name);
  if (e) {
    return atof(e);
  } else {
    return def;
  }
}

static int envr(mpfr_t out, const char *name, const char *def) {
  const char *e = getenv(name);
  if (e) {
    return mpfr_set_str(out, e,   10, MPFR_RNDN);
  } else {
    return mpfr_set_str(out, def, 10, MPFR_RNDN);
  }
}

extern int main(int argc, char **argv) {
  int workers = envi("threads", 4);
  int width = envi("width", 1280);
  int height = envi("height", 720);
  int detect_glitches = envi("detect_glitches", 1);
  int approx_skip = envi("approx_skip", 0);
  int maxiters = envi("maxiters", 1 << 18);
  int chunk = envi("chunk", 1 << 8);
  double escape_radius = envd("escaperadius", 25);
  double glitch_threshold = envd("glitchthreshold", 1e-6);
  int precision = envi("precision", 53);
  mpfr_t radius, centerx, centery;
  mpfr_init2(radius, 53);
  envr(radius, "radius", "2.0");
  int e = max(53, 53 - mpfr_get_exp(radius));
  if (e > precision) {
    fprintf(stderr, "WARNING: increasing precision to %d\n", e);
    precision = e;
  }
  mpfr_init2(centerx, precision);
  mpfr_init2(centery, precision);
  envr(centerx, "real", "-0.75");
  envr(centery, "imag", "0.0");
  struct perturbator *context = perturbator_new(workers, width, height, maxiters, chunk, escape_radius, glitch_threshold);
  perturbator_set_detect_glitches(context, detect_glitches);
  perturbator_set_approx_skip(context, approx_skip);
  perturbator_start(context, centerx, centery, radius);
  perturbator_stop(context, false);
  const float *data = perturbator_get_output(context);
  fwrite(data, 4 * width * height * sizeof(float), 1, stdout);
  fflush(stdout);
  return 0;
(void) argc;
(void) argv;
}
