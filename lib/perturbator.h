// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015,2016 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#ifndef PERTURBATOR_H
#define PERTURBATOR_H 1

#include <mpfr.h>

#ifdef __cplusplus
extern "C" {
#endif

struct perturbator;
struct perturbator *perturbator_new(int workers, int width, int height, int maxiters, int chunk, double escape_radius, double glitch_threshold);
void perturbator_start(struct perturbator *context, const mpfr_t x, const mpfr_t y, const mpfr_t r);
void perturbator_stop(struct perturbator *context, int force);
int perturbator_active(struct perturbator *context);
const float *perturbator_get_output(struct perturbator *context);
int perturbator_get_primary_reference(struct perturbator *context, mpfr_t x, mpfr_t y);
void perturbator_set_detect_glitches(struct perturbator *img, int detect_glitches);
void perturbator_set_approx_skip(struct perturbator *img, int approx_skip);

#ifdef __cplusplus
}
#endif

#endif
