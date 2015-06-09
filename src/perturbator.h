#ifndef PERTURBATOR_H
#define PERTURBATOR_H 1

#include <mpfr.h>

struct perturbator;
struct perturbator *perturbator_new(int workers, int width, int height, int maxiters, int chunk, double escape_radius, double glitch_threshold, int precision);
void perturbator_start(struct perturbator *context, const mpfr_t x, const mpfr_t y, const mpfr_t r);
void perturbator_stop(struct perturbator *context, int force);
int perturbator_active(struct perturbator *context);
const float *perturbator_get_output(struct perturbator *context);

#endif
