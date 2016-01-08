// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015,2016 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#include <assert.h>
#include <complex>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>

#include "perturbator.h"

#include "edouble.cc"

#include "generated/z2c.h"

extern "C" {
  using namespace std;
#include "mandelbrot-numerics.h"
#include "m_r_nucleus.c"
#include "m_r_shape.c"
#include "m_r_box_period.c"
#include "m_r_domain_size.c"
}

enum float_type
  { ft_float
  , ft_double
  , ft_long_double
  , ft_edouble
  };

// don't use float, derivative overflows too easily...
//#define EXP_THRESHOLD_FLOAT -120
#define EXP_THRESHOLD_DOUBLE -960
#define EXP_THRESHOLD_LONG_DOUBLE -16300
#define EXP_THRESHOLD_EDOUBLE -(1<<29)

struct series_node {
  struct series_node *next;
  int exponent;
  int iters;
  mpc_t z;
  enum float_type ft;
  union {
    struct z2c_approx<float> *approxf;
    struct z2c_approx<double>  *approx;
    struct z2c_approx<long double> *approxl;
    struct z2c_approx<edouble> *approxe;
  } u;
};


template <typename R>
struct reference;

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
  void *volatile refs;

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


struct series_node *image_cached_approx(struct perturbator *img, bool reused, const mpc_t c, int exponent, mpc_t z, int *iter);


float mpfr_get(const mpfr_t &op, mpfr_rnd_t rnd, float dummy) {
  (void) dummy;
  return mpfr_get_flt(op, rnd);
}

double mpfr_get(const mpfr_t &op, mpfr_rnd_t rnd, double dummy) {
  (void) dummy;
  return mpfr_get_d(op, rnd);
}

long double mpfr_get(const mpfr_t &op, mpfr_rnd_t rnd, long double dummy) {
  (void) dummy;
  return mpfr_get_ld(op, rnd);
}

edouble mpfr_get(const mpfr_t &op, mpfr_rnd_t rnd, edouble dummy) {
  (void) dummy;
  (void) rnd;
  return edouble(op);
}

int mpfr_add(mpfr_t &rop, const mpfr_t &op1, float op2, mpfr_rnd_t rnd) {
  return mpfr_add_d(rop, op1, op2, rnd);
}

int mpfr_add(mpfr_t &rop, const mpfr_t &op1, double op2, mpfr_rnd_t rnd) {
  return mpfr_add_d(rop, op1, op2, rnd);
}

int mpfr_add(mpfr_t &rop, const mpfr_t &op1, long double op2, mpfr_rnd_t rnd) {
  mpfr_t tmp;
  mpfr_init2(tmp, 64);
  to_mpfr(op2, tmp);
  int r = mpfr_add(rop, op1, tmp, rnd);
  mpfr_clear(tmp);
  return r;
}

int mpfr_add(mpfr_t &rop, const mpfr_t &op1, edouble op2, mpfr_rnd_t rnd) {
  mpfr_t tmp;
  mpfr_init2(tmp, 53);
  to_mpfr(op2, tmp);
  int r = mpfr_add(rop, op1, tmp, rnd);
  mpfr_clear(tmp);
  return r;
}

template <typename T>
void to_mpc(const std::complex<T> &from, mpc_t &to) {
  to_mpfr(std::real(from), mpc_realref(to));
  to_mpfr(std::imag(from), mpc_imagref(to));
}

template <typename T>
int mpc_add(mpc_t &rop, const mpc_t &op1, const std::complex<T> &op2, mpc_rnd_t rnd) {
  (void) rnd;
          mpfr_add(mpc_realref(rop), mpc_realref(op1), std::real(op2), MPFR_RNDN);
  return  mpfr_add(mpc_imagref(rop), mpc_imagref(op1), std::imag(op2), MPFR_RNDN);
}

template <typename T>
int mpc_div(mpc_t &rop, const mpc_t &op1, const std::complex<T> &op2, mpc_rnd_t rnd) {
  mpc_t tmp;
  mpc_init2(tmp, 64);
  to_mpc(op2, tmp);
  int r = mpc_div(rop, op1, tmp, rnd);
  mpc_clear(tmp);
  return r;
}

template <typename T>
std::complex<T> mpc_get(const mpc_t &op, mpc_rnd_t rnd, T dummy) {
  (void) rnd;
  return std::complex<T>(mpfr_get(mpc_realref(op), MPFR_RNDN, dummy), mpfr_get(mpc_imagref(op), MPFR_RNDN, dummy));
}

#include "core_native.c"

extern struct perturbator *perturbator_new(int workers, int width, int height, int maxiters, int chunk, double escape_radius, double glitch_threshold) {
  struct perturbator *img = (struct perturbator *) calloc(1, sizeof(*img));

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

  img->output = (float *) calloc(1, width * height * 4 * sizeof(*img->output));

  pthread_mutex_init(&img->mutex, 0);
  pthread_cond_init(&img->cond, 0);

  image_log(img, LOG_QUEUE | LOG_CACHE, "sequence action   queued   iters  period    active count\n");
  image_log(img, LOG_QUEUE | LOG_CACHE, "-------- ------ -------- ------- ------- ---------------\n");
  return img;
}


extern void perturbator_start(struct perturbator *img, const mpfr_t centerx, const mpfr_t centery, const mpfr_t radius) {

  img->precision = max(53, int(53 - 2 * mpfr_get_exp(radius)));

  image_log(img, LOG_VIEW, "real=%Re\nimag=%Re\nradius=%Re\nprecision=%d\n", centerx, centery, radius, img->precision);

  mpc_set_prec(img->center, img->precision);
  mpc_set_fr_fr(img->center, centerx, centery, MPC_RNDNN);
  mpfr_set(img->radius, radius, MPFR_RNDN);

  memset(img->output, 0, img->width * img->height * 4 * sizeof(*img->output));

  int e = mpfr_get_exp(radius);
  img->ft =
    e >= EXP_THRESHOLD_DOUBLE ? ft_double :
    e >= EXP_THRESHOLD_LONG_DOUBLE ? ft_long_double :
    ft_edouble;
  switch (img->ft) {
    case ft_float:        perturbator_start_internal<float>(img); return;
    case ft_double:       perturbator_start_internal<double>(img); return;
    case ft_long_double:  perturbator_start_internal<long double>(img); return;
    case ft_edouble:      perturbator_start_internal<edouble>(img); return;
  }
  assert(! "valid float type");
}


extern void perturbator_stop(struct perturbator *img, int force) {
  if (force) {
    img->running = false;
  }
  pthread_mutex_lock(&img->mutex);
  pthread_cond_broadcast(&img->cond);
  pthread_mutex_unlock(&img->mutex);
  for (int i = 0; i < img->workers; ++i) {
    pthread_join(img->threads[i], 0);
  }
  switch (img->ft) {
    case ft_float:        release_refs<float>(img); return;
    case ft_double:       release_refs<double>(img); return;
    case ft_long_double:  release_refs<long double>(img); return;
    case ft_edouble:      release_refs<edouble>(img); return;
  }
  assert(! "valid float type");
}


void image_delete_cache(struct perturbator *img) {
  if (img->series) {
    z2c_series_delete(img->series);
    img->series = 0;
  }
  while (img->nodes) {
    struct series_node *next = img->nodes->next;
    mpc_clear(img->nodes->z);
    switch (img->nodes->ft) {
      case ft_float:       z2c_approx_delete(img->nodes->u.approxf); break;
      case ft_double:      z2c_approx_delete(img->nodes->u.approx ); break;
      case ft_long_double: z2c_approx_delete(img->nodes->u.approxl); break;
      case ft_edouble:     z2c_approx_delete(img->nodes->u.approxe); break;
      default: assert(! "valid float type");
    }
    img->nodes = next;
  }
}

struct series_node *image_cached_approx(struct perturbator *img, bool reused, const mpc_t c, int exponent, mpc_t z, int *iter) {
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
        struct series_node *node = (struct series_node *) calloc(1, sizeof(*node));
        node->next = img->nodes;
        node->exponent = e;
        node->iters = iters;
        mpc_init2(node->z, 53); // updated by get_z
        struct z2c_reference *reference = z2c_series_reference_new(img->series);
        z2c_reference_get_zr(reference, mpc_realref(node->z), mpc_imagref(node->z));
        z2c_reference_delete(reference);
        // don't use float, derivative overflows too easily...
        /*if (e >= EXP_THRESHOLD_FLOAT) {
          node->ft = ft_float;
          node->u.approxf = z2c_series_approx_newf(img->series, e);
        } else*/ if (e >= EXP_THRESHOLD_DOUBLE) {
          node->ft = ft_double;
          node->u.approx  = z2c_series_approx_new(img->series, e, double(0));
        } else if (e >= EXP_THRESHOLD_LONG_DOUBLE) {
          node->ft = ft_long_double;
          node->u.approxl = z2c_series_approx_new(img->series, e, (long double)(0));
        } else if (e >= EXP_THRESHOLD_EDOUBLE) {
          node->ft = ft_edouble;
          node->u.approxe = z2c_series_approx_new(img->series, e, edouble(0));
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
