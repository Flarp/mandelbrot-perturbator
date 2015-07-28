// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>

#include "perturbator.h"
#include "generated/z2c.h"

#include "mandelbrot-numerics.h"
#include "m_r_nucleus.c"
#include "m_r_shape.c"
#include "m_r_box_period.c"
#include "m_r_domain_size.c"

static inline int max(int a, int b) {
  return a > b ? a : b;
}

enum float_type
  { ft_float
  , ft_double
  , ft_long_double
  };

//#define EXP_THRESHOLD_FLOAT -120
#define EXP_THRESHOLD_DOUBLE -1016
#define EXP_THRESHOLD_LONG_DOUBLE -16376

struct series_node {
  struct series_node *next;
  int exponent;
  int iters;
  mpc_t z;
  enum float_type ft;
  union {
    struct z2c_approxf *approxf;
    struct z2c_approx  *approx;
    struct z2c_approxl *approxl;
  } u;
};


struct referencef;
struct reference;
struct referencel;

struct perturbator {
  int workers;
  int width;
  int height;
  int maxiters;
  double escape_radius;
  double glitch_threshold;
  int precision;
  int chunk;
  mpc_t center;
  mpfr_t radius;
  double escape_radius_2;
  double log_escape_radius_2;
  int newton_steps_root;
  int newton_steps_child;
  int order;
  int threshold;
  int logging;

  // cache
  mpc_t last_reference;
  int last_period;
  struct z2c_series *series;
  struct series_node *nodes;

  float *output;

  pthread_mutex_t mutex;
  pthread_cond_t cond;
  pthread_t *threads;

  int volatile active_workers;
  bool volatile running;

  int volatile queue_id;
  int volatile start_id;

  enum float_type ft;
  union {
    struct referencef *volatile refsf;
    struct reference  *volatile refs;
    struct referencel *volatile refsl;
  } urefs;

};

extern const float *perturbator_get_output(struct perturbator *img) {
  return img->output;
}


static bool image_running(struct perturbator *img) {
  return img->running;
}

extern int perturbator_active(struct perturbator *img) {
  pthread_mutex_lock(&img->mutex);
  bool active = img->active_workers;
  pthread_mutex_unlock(&img->mutex);
  return active;
}

#define LOG_VIEW  1
#define LOG_QUEUE 2
#define LOG_CACHE 4

#define image_log(img, aspect, ...) do{ \
  if ((img)->logging & (aspect)) { mpfr_fprintf(stderr, __VA_ARGS__); } \
  }while(0)


static struct series_node *image_cached_approx(struct perturbator *img, bool reused, const mpc_t c, int exponent, mpc_t z, int *iter);

// FIXME
static void mpfr_add_ld(mpfr_t rop, const mpfr_t op1, long double op2, mpfr_rnd_t rnd) {
  mpfr_t tmp;
  mpfr_init2(tmp, 64);
  mpfr_set_ld(tmp, op2, rnd);
  mpfr_add(rop, op1, tmp, rnd);
  mpfr_clear(tmp);
}

#define FT ft_float
#define FTYPE float
#define FNAME(x) x ## f
#define FMPFRGET mpfr_get_flt
#define FMPFRADD mpfr_add_d
#define FMPCGET mpc_get_fltc
#define FMPCADD mpc_add_fltc
#define FMPCDIV mpc_div_fltc2
#include "core_native.c"
#undef FTYPE
#undef FNAME
#undef FMPFRGET
#undef FMPFRADD
#undef FMPCGET
#undef FMPCADD
#undef FMPCDIV
#undef FT

#define FT ft_double
#define FTYPE double
#define FNAME(x) x
#define FMPFRGET mpfr_get_d
#define FMPFRADD mpfr_add_d
#define FMPCGET mpc_get_dc2
#define FMPCADD mpc_add_dc
#define FMPCDIV mpc_div_dc2
#include "core_native.c"
#undef FTYPE
#undef FNAME
#undef FMPFRGET
#undef FMPFRADD
#undef FMPCGET
#undef FMPCADD
#undef FMPCDIV
#undef FT

#define FT ft_long_double
#define FTYPE long double
#define FNAME(x) x ## l
#define FMPFRGET mpfr_get_ld
#define FMPFRADD mpfr_add_ld
#define FMPCGET mpc_get_ldc2
#define FMPCADD mpc_add_ldc2
#define FMPCDIV mpc_div_ldc2
#include "core_native.c"
#undef FTYPE
#undef FNAME
#undef FMPFRGET
#undef FMPFRADD
#undef FMPCGET
#undef FMPCADD
#undef FMPCDIV
#undef FT

extern struct perturbator *perturbator_new(int workers, int width, int height, int maxiters, int chunk, double escape_radius, double glitch_threshold) {
  struct perturbator *img = calloc(1, sizeof(*img));

  img->workers = workers;
  img->width = width;
  img->height = height;
  img->maxiters = maxiters;
  img->chunk = chunk;
  img->escape_radius = escape_radius;
  img->glitch_threshold = glitch_threshold;
  img->precision = 53;

  mpc_init2(img->center, 53);
  mpfr_init2(img->radius, 53);

  img->escape_radius_2 = escape_radius * escape_radius;
  img->log_escape_radius_2 = log(img->escape_radius_2);

  img->newton_steps_root = 64;
  img->newton_steps_child = 8;

  img->order = 24;
  img->threshold = 64;
  img->logging = -1;

  mpc_init2(img->last_reference, 53);
  mpc_set_ui_ui(img->last_reference, 0, 0, MPC_RNDNN);
  img->last_period = 1;

  img->output = calloc(1, width * height * 4 * sizeof(*img->output));

  pthread_mutex_init(&img->mutex, 0);
  pthread_cond_init(&img->cond, 0);

  image_log(img, LOG_QUEUE | LOG_CACHE, "sequence action   queued   iters  period    active count\n");
  image_log(img, LOG_QUEUE | LOG_CACHE, "-------- ------ -------- ------- ------- ---------------\n");
  return img;
}


extern void perturbator_start(struct perturbator *img, const mpfr_t centerx, const mpfr_t centery, const mpfr_t radius) {

  img->precision = max(53, 53 - 2 * mpfr_get_exp(radius));

  image_log(img, LOG_VIEW, "real=%Re\nimag=%Re\nradius=%Re\nprecision=%d\n", centerx, centery, radius, img->precision);

  mpc_set_prec(img->center, img->precision);
  mpc_set_fr_fr(img->center, centerx, centery, MPC_RNDNN);
  mpfr_set(img->radius, radius, MPFR_RNDN);

  memset(img->output, 0, img->width * img->height * 4 * sizeof(*img->output));

  img->ft = mpfr_get_exp(radius) >= EXP_THRESHOLD_DOUBLE ? ft_double : ft_long_double;
  switch (img->ft) {
    case ft_float:        perturbator_start_internalf(img); return;
    case ft_double:       perturbator_start_internal (img); return;
    case ft_long_double:  perturbator_start_internall(img); return;
  }
  assert(! "valid float type");
}


extern void perturbator_stop(struct perturbator *img, int force) {
  pthread_mutex_lock(&img->mutex);
  if (force) {
    img->running = false;
  }
  pthread_cond_broadcast(&img->cond);
  pthread_mutex_unlock(&img->mutex);
  for (int i = 0; i < img->workers; ++i) {
    pthread_join(img->threads[i], 0);
  }
  switch (img->ft) {
    case ft_float:        release_refsf(img); return;
    case ft_double:       release_refs (img); return;
    case ft_long_double:  release_refsl(img); return;
  }
  assert(! "valid float type");
}


static void image_delete_cache(struct perturbator *img) {
  if (img->series) {
    z2c_series_delete(img->series);
    img->series = 0;
  }
  while (img->nodes) {
    struct series_node *next = img->nodes->next;
    mpc_clear(img->nodes->z);
    switch (img->nodes->ft) {
      case ft_float:       z2c_approx_deletef(img->nodes->u.approxf); break;
      case ft_double:      z2c_approx_delete (img->nodes->u.approx ); break;
      case ft_long_double: z2c_approx_deletel(img->nodes->u.approxl); break;
      default: assert(! "valid float type");
    }
    img->nodes = next;
  }
}

static struct series_node *image_cached_approx(struct perturbator *img, bool reused, const mpc_t c, int exponent, mpc_t z, int *iter) {
  pthread_mutex_lock(&img->mutex);
  if (! reused) {
    image_delete_cache(img);
  }
  if (! img->series) {
    img->series = z2c_series_new(img->order, mpc_realref(c), mpc_imagref(c));
  }
  int exponent0 = 16;
  if (img->nodes)  {
    exponent0 = img->nodes->exponent;
  }
  // cache all the things
  for (int e = exponent0; e >= exponent; --e) {
    while (z2c_series_step(img->series, e + 2, img->threshold)) {
//      if (! image_running(img)) {
//        return 0;
//      }
      // nop
    }
    int iters = z2c_series_get_n(img->series);
    if ((! img->nodes) || (img->nodes && iters >= img->nodes->iters)) {
      if (img->nodes && iters == img->nodes->iters) {
        img->nodes->exponent = e;
      } else {
        struct series_node *node = calloc(1, sizeof(*node));
        node->next = img->nodes;
        node->exponent = e;
        node->iters = iters;
        mpc_init2(node->z, 53); // updated by get_z
        struct z2c_reference *reference = z2c_series_reference_new(img->series);
        z2c_reference_get_zr(reference, mpc_realref(node->z), mpc_imagref(node->z));
        z2c_reference_delete(reference);
        /*if (e >= EXP_THRESHOLD_FLOAT) {
          node->ft = ft_float;
          node->u.approxf = z2c_series_approx_newf(img->series, e);
        } else*/ if (e >= EXP_THRESHOLD_DOUBLE) {
          node->ft = ft_double;
          node->u.approx  = z2c_series_approx_new (img->series, e);
        } else if (e >= EXP_THRESHOLD_LONG_DOUBLE) {
          node->ft = ft_long_double;
          node->u.approxl = z2c_series_approx_newl(img->series, e);
        } else {
          assert(! "exponent in range of supported types");
        }
        img->nodes = node;
      }
    }
  }
  // lookup in the cache
  struct series_node *node = img->nodes, *prev = img->nodes;
  while (node && node->exponent < exponent) {
    prev = node;
    node = node->next;
  }
  if (node) {
    mpc_set(z, node->z, MPC_RNDNN);
    *iter = node->iters;
    pthread_mutex_unlock(&img->mutex);
    return node;
  } else {
    pthread_mutex_unlock(&img->mutex);
    return prev;
  }
}

int perturbator_get_primary_reference(struct perturbator *img, mpfr_t x, mpfr_t y) {
  pthread_mutex_lock(&img->mutex);
  int period = img->last_period;
  int prec = mpc_get_prec(img->last_reference);
  mpfr_set_prec(x, prec);
  mpfr_set_prec(y, prec);
  mpfr_set(x, mpc_realref(img->last_reference), MPFR_RNDN);
  mpfr_set(y, mpc_imagref(img->last_reference), MPFR_RNDN);
  pthread_mutex_unlock(&img->mutex);
  return period;
}
