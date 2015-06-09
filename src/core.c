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


static inline double cnorm(complex double z) {
  return creal(z) * creal(z) + cimag(z) * cimag(z);
}


static inline void mpc_add_dc(mpc_t rop, const mpc_t op1, complex double op2, mpc_rnd_t rnd) {
(void) rnd;
  mpfr_add_d(mpc_realref(rop), mpc_realref(op1), creal(op2), MPFR_RNDN);
  mpfr_add_d(mpc_imagref(rop), mpc_imagref(op1), cimag(op2), MPFR_RNDN);
}


struct pixel {
  complex double c;
  complex double z;
  uint32_t index;
  uint32_t iters;
};


static int cmp_pixel_by_iters_asc(const void *a, const void *b) {
  const struct pixel *x = a;
  const struct pixel *y = b;
  if (x->iters < y->iters) { return -1; }
  if (x->iters > y->iters) { return  1; }
  return 0;
}


struct reference {
  struct reference *next;
  struct reference *parent;
  int use_count;
  int queue_id;
  int start_id;

  int period;
  int iters;
  mpc_t c;
  mpc_t z;
  complex double z_d_old;
  complex double z_d;
  double z_d2eM6;

  struct pixel *px[2];
  int count;

  int index;
};


struct series_node {
  struct series_node *next;
  int exponent;
  int iters;
  mpc_t z;
  struct z2c_approx *approx;
};


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
  struct reference *volatile refs;

};


extern const float *perturbator_get_output(struct perturbator *img) {
  return img->output;
}


static bool image_running(struct perturbator *img) {
  pthread_mutex_lock(&img->mutex);
  bool running = img->running;
  pthread_mutex_unlock(&img->mutex);
  return running;
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


static void reference_release(struct reference *ref) {
  if (ref) {
    ref->use_count -= 1;
    if (ref->use_count <= 0) {
      reference_release(ref->parent);
      ref->parent = 0;
      if (ref->px[0]) {
        free(ref->px[0]);
        ref->px[0] = 0;
      }
      if (ref->px[1]) {
        free(ref->px[1]);
        ref->px[1] = 0;
      }
      mpc_clear(ref->c);
      mpc_clear(ref->z);
      free(ref);
    }
  }
}


static struct reference *image_dequeue(struct perturbator *img, struct reference *ref) {
  pthread_mutex_lock(&img->mutex);
  if (ref) {
    image_log(img, LOG_QUEUE, "%8d DONE   %8d\n", ref->start_id, ref->queue_id);
  }
  // release the old reference
  reference_release(ref);
  // if no workers are working and the queue is empty, it would be empty forever
  img->active_workers -= 1;
  while (img->running && ! img->refs) {
    if (! img->active_workers) {
      img->running = false;
      pthread_cond_broadcast(&img->cond);
    } else {
      pthread_cond_wait(&img->cond, &img->mutex);
    }
  }

  if (img->running && img->refs) {
    // still work to do
    ref = img->refs;
    img->refs = ref->next;
    ref->start_id = (img->start_id += 1);
    ref->use_count += 1;
    img->active_workers += 1;
    image_log(img, LOG_QUEUE, "%8d START  %8d%8d%8d%16d\n", ref->start_id, ref->queue_id, ref->iters, ref->period, ref->count);
  } else {
    ref = 0;
  }
  pthread_mutex_unlock(&img->mutex);
  return ref;
}


static void image_enqueue(struct perturbator *img, struct reference *ref) {
  pthread_mutex_lock(&img->mutex);
  ref->queue_id = (img->queue_id += 1);
  image_log(img, LOG_QUEUE, "         QUEUE  %8d%8d%8d%16d\n", ref->queue_id, ref->iters, ref->period, ref->count);
  // ensure parent isn't deallocated
  if (ref->parent) {
    ref->parent->use_count += 1;
    mpc_init2(ref->z, img->precision);
    mpc_init2(ref->c, img->precision);
  }
  // link into reference queue (sorted by count descending)
  struct reference *before = 0;
  struct reference *after = img->refs;
  while (after && after->count > ref->count) {
    before = after;
    after = after->next;
  }
  ref->next = after;
  if (before) {
    before->next = ref;
  } else {
    img->refs = ref;
  }
  // wake waiting workers
  pthread_cond_broadcast(&img->cond);
  pthread_mutex_unlock(&img->mutex);
}


static void *image_worker(void *arg);


extern struct perturbator *perturbator_new(int workers, int width, int height, int maxiters, int chunk, double escape_radius, double glitch_threshold, int precision) {
  struct perturbator *img = calloc(1, sizeof(*img));

  img->workers = workers;
  img->width = width;
  img->height = height;
  img->maxiters = maxiters;
  img->chunk = chunk;
  img->escape_radius = escape_radius;
  img->glitch_threshold = glitch_threshold;
  img->precision = precision;
  
  mpc_init2(img->center, precision);
  mpfr_init2(img->radius, 53);

  img->escape_radius_2 = escape_radius * escape_radius;
  img->log_escape_radius_2 = log(img->escape_radius_2);

  img->newton_steps_root = 64;
  img->newton_steps_child = 8;

  img->order = 24;
  img->threshold = 64;
  img->logging = LOG_VIEW;

  mpc_init2(img->last_reference, precision);
  mpc_set_ui_ui(img->last_reference, 0, 0, MPC_RNDNN);

  img->output = calloc(1, width * height * 4 * sizeof(*img->output));

  pthread_mutex_init(&img->mutex, 0);
  pthread_cond_init(&img->cond, 0);

  image_log(img, LOG_QUEUE | LOG_CACHE, "sequence action   queued   iters  period    active count\n");
  image_log(img, LOG_QUEUE | LOG_CACHE, "-------- ------ -------- ------- ------- ---------------\n");
  return img;
}


extern void perturbator_start(struct perturbator *img, const mpfr_t centerx, const mpfr_t centery, const mpfr_t radius) {

  image_log(img, LOG_VIEW, "real=%Re\nimag=%Re\nradius=%Re\n", centerx, centery, radius);

  mpc_set_fr_fr(img->center, centerx, centery, MPC_RNDNN);
  mpfr_set(img->radius, radius, MPFR_RNDN);

  memset(img->output, 0, img->width * img->height * 4 * sizeof(*img->output));

  struct reference *ref = calloc(1, sizeof(*ref));
  mpc_init2(ref->c, img->precision);
  mpc_init2(ref->z, img->precision);
  mpc_set(ref->c, img->center, MPC_RNDNN);
  mpc_set_ui_ui(ref->z, 0, 0, MPC_RNDNN);
  ref->iters = 0;
  ref->period = 0;
  ref->count = img->width * img->height;
  ref->px[0] = malloc(ref->count * sizeof(*ref->px[0]));
  double vdiameter = 2.0 * mpfr_get_d(img->radius, MPFR_RNDN);
  double hdiameter = img->width * vdiameter / img->height;
  #pragma omp parallel for
  for (int j = 0; j < img->height; ++j) {
    double y = ((j + 0.5) / img->height - 0.5) * vdiameter;
    for (int i = 0; i < img->width; ++i) {
      double x = ((i + 0.5) / img->width - 0.5) * hdiameter;
      complex double c = x - I * y;
      int index = j * img->width + i;
      ref->px[0][index].c = c;
      ref->px[0][index].z = 0;
      ref->px[0][index].index = index;
    }
  }
  img->threads = calloc(1, img->workers * sizeof(*img->threads));
  img->running = true;
  image_enqueue(img, ref);
  pthread_mutex_lock(&img->mutex);
  img->active_workers = img->workers;
  for (int i = 0; i < img->workers; ++i) {
    pthread_create(&img->threads[i], 0, image_worker, img);
  }
  pthread_mutex_unlock(&img->mutex);
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
  for (struct reference *ref = img->refs; ref; ) {
    struct reference *next = ref->next;
    reference_release(ref);
    ref = next;
  }
  img->refs = 0;
}


static void image_delete_cache(struct perturbator *img) {
  if (img->series) {
    z2c_series_delete(img->series);
    img->series = 0;
  }
  while (img->nodes) {
    struct series_node *next = img->nodes->next;
    mpc_clear(img->nodes->z);
    z2c_approx_delete(img->nodes->approx);
    img->nodes = next;
  }
}


static struct z2c_approx *image_cached_approx(struct perturbator *img, bool reused, const mpc_t c, int exponent, mpc_t z, int *iter) {
  pthread_mutex_lock(&img->mutex);
  if (! reused) {
    image_delete_cache(img);
  }
  if (! img->series) {
    img->series = z2c_series_new(img->order, mpc_realref(c), mpc_imagref(c));
  }
  int exponent0 = 8;
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
        mpc_init2(node->z, img->precision);
        struct z2c_reference *reference = z2c_series_reference_new(img->series);
        z2c_reference_get_zr(reference, mpc_realref(node->z), mpc_imagref(node->z));
        z2c_reference_delete(reference);
        node->approx = z2c_series_approx_new(img->series, e);
        img->nodes = node;
      }
    }
  }
  // lookup in the cache
  struct series_node *node = img->nodes;
  while (node && node->exponent < exponent) {
    node = node->next;
  }
  if (node) {
    mpc_set(z, node->z, MPC_RNDNN);
    *iter = node->iters;
    pthread_mutex_unlock(&img->mutex);
    return node->approx;
  } else {
    pthread_mutex_unlock(&img->mutex);
    return 0;
  }
}


static void *image_worker(void *arg) {
  struct perturbator *img = arg;
  image_log(img, LOG_QUEUE, "         ENTER\n");
  int maxiters = img->maxiters;
  double escape_radius_2 = img->escape_radius_2;
  double log_escape_radius_2 = img->log_escape_radius_2;
  int newton_steps_root = img->newton_steps_root;
  int newton_steps_child = img->newton_steps_child;
  int chunk = img->chunk;
  double glitch_threshold = img->glitch_threshold;
  int precision = img->precision;

  complex double *z_d = malloc((1 + chunk) * sizeof(*z_d));
  double *z_size = malloc((1 + chunk) * sizeof(*z_size));

  struct reference *ref = 0;
  while ( (ref = image_dequeue(img, ref)) ) {

    // find reference atom
    if (ref->parent) {

      mpc_add_dc(ref->c, ref->parent->c, ref->px[0][ref->index].c, MPC_RNDNN);
      m_r_nucleus(ref->c, ref->c, ref->period, newton_steps_child);
      mpc_sub(ref->z, ref->parent->c, ref->c, MPC_RNDNN);
      complex double dc = -mpc_get_dc(ref->z, MPC_RNDNN);
      mpc_set_ui_ui(ref->z, 0, 0, MPC_RNDNN);
      // rebase pixels to new reference
      int count = ref->count;
      for (int k = 0; k < count; ++k) {
        ref->px[0][k].c += dc;
      }

    } else {

      int exponent = mpfr_get_exp(img->radius);
      // find an appropriate initial reference
      bool ok = false;
      bool reused = false;
      mpc_t nucleus;
      mpc_init2(nucleus, precision);
      mpc_set_ui_ui(nucleus, 0, 0, MPC_RNDNN);
      if (! ok) {
        // try reuse
        mpc_t delta;
        mpc_init2(delta, 53);
        mpfr_t delta2, radius2;
        mpfr_init2(delta2, 53);
        mpfr_init2(radius2, 53);
        mpfr_sqr(radius2, img->radius, MPFR_RNDN);
        mpfr_mul_d(radius2, radius2, 65536, MPFR_RNDN);
        mpc_sub(delta, img->last_reference, img->center, MPC_RNDNN);
        mpc_norm(delta2, delta, MPFR_RNDN);
        if (mpfr_less_p(delta2, radius2)) {
          mpc_set(nucleus, img->last_reference, MPC_RNDNN);
          if (! mpfr_zero_p(delta2)) {
            exponent = fmax(exponent, mpfr_get_exp(delta2) / 2);
          }
          ok = true;
          reused = true;
          image_log(img, LOG_CACHE, "         REUSE\n");
        }
        mpc_clear(delta);
        mpfr_clear(delta2);
        mpfr_clear(radius2);
      }
      if (! ok) {
        // try box period
        int period = m_r_box_period_do(img->center, img->radius, maxiters);
        if (period) {
          m_r_nucleus(nucleus, img->center, period, newton_steps_root);
          if (m_cardioid == m_r_shape(nucleus, period)) {
            ok = true;
            image_log(img, LOG_CACHE, "         BOXED                  %8d\n", period);
          }
        }
      }
      if (! ok) {
        // try modified partials
        image_log(img, LOG_CACHE, "         SEARCH\n");
        mpc_t z, dc0;
        mpc_init2(z, precision);
        mpc_set_ui_ui(z, 0, 0, MPC_RNDNN);
        mpc_init2(dc0, 53);
        mpfr_t z2, mz2, dc2, radius2;
        mpfr_init2(z2, 53);
        mpfr_init2(mz2, 53);
        mpfr_set_d(mz2, 65536, MPFR_RNDN);
        mpfr_init2(dc2, 53);
        mpfr_init2(radius2, 53);
        mpfr_sqr(radius2, img->radius, MPFR_RNDN);
        mpfr_mul_d(radius2, radius2, 65536, MPFR_RNDN);

        for (int period = 1; period < maxiters; ++period) {
          if (! image_running(img)) {
            break;
          }
          mpc_sqr(z, z, MPC_RNDNN);
          mpc_add(z, z, img->center, MPC_RNDNN);
          mpc_norm(z2, z, MPFR_RNDN);
          if (mpfr_get_d(z2, MPFR_RNDN) > escape_radius_2) {
            break;
          }
          mpfr_mul_2si(z2, z2, -period, MPFR_RNDN);
          if (mpfr_less_p(z2, mz2)) {
            mpfr_set(mz2, z2, MPFR_RNDN);
            m_r_nucleus(ref->z, ref->c, period, newton_steps_root);
            if (m_cardioid == m_r_shape(ref->z, period)) {
              mpc_sub(dc0, ref->c, ref->z, MPC_RNDNN);
              mpc_norm(dc2, dc0, MPFR_RNDN);
              if (mpfr_less_p(dc2, radius2)) {
                mpc_set(nucleus, ref->z, MPC_RNDNN);
                ok = true;
                image_log(img, LOG_CACHE, "         FOUND                  %8d\n", period);
                break;
              }
            }
          }
        }
        mpc_clear(z);
        mpc_clear(dc0);
        mpfr_clear(z2);
        mpfr_clear(mz2);
        mpfr_clear(dc2);
        mpfr_clear(radius2);
      }
      if (! ok && image_running(img)) {
        image_log(img, LOG_CACHE, "         REUSE FALLBACK\n");
      }

      if (ok) {
        mpc_set(img->last_reference, nucleus, MPC_RNDNN);
      }
      mpc_set(ref->c, nucleus, MPC_RNDNN);
      mpc_sub(nucleus, ref->c, img->center, MPC_RNDNN);
      complex double dc = -mpc_get_dc(nucleus, MPC_RNDNN);
      mpc_clear(nucleus);
      mpc_set_ui_ui(ref->z, 0, 0, MPC_RNDNN);
      ref->z_d_old = 0;
      ref->z_d = 0;
      struct z2c_approx *approx = image_cached_approx(img, reused, ref->c, exponent, ref->z, &ref->iters);
      ref->z_d = mpc_get_dc(ref->z, MPC_RNDNN);

      // rebase pixels to new reference with approximation
      int count = ref->count;
      #pragma omp parallel for
      for (int k = 0; k < count; ++k) {
        ref->px[0][k].c += dc;
        ref->px[0][k].z = z2c_approx_do(approx, ref->px[0][k].c);
      }

      image_log(img, LOG_CACHE, "         APPROX         %8d\n", ref->iters);

    } // if parent

    // prepare for output pixels
    ref->px[1] = malloc(ref->count * sizeof(*ref->px[1]));
    struct pixel *glitched = malloc(ref->count * sizeof(*glitched));
    int start_id = ref->start_id;

    for (int iters = ref->iters; iters < maxiters; iters += chunk) {
      if (! image_running(img)) {
        break;
      }

      // compute ref chunk
      z_d[0] = ref->z_d;
      for (int k = 1; k <= chunk; ++k) {
        mpc_sqr(ref->z, ref->z, MPC_RNDNN);
        mpc_add(ref->z, ref->z, ref->c, MPC_RNDNN);
        z_d[k] = mpc_get_dc(ref->z, MPC_RNDNN);
      }
      ref->z_d = z_d[chunk];
      for (int k = 0; k <= chunk; ++k) {
        z_size[k] = cnorm(z_d[k]) * glitch_threshold;
      }

      // advance pixels
      int input_count = ref->count;
      int active_count = 0;
      int glitch_count = 0;

      #pragma omp parallel for      
      for (int k = 0; k < input_count; ++k) {
        struct pixel *in = &ref->px[0][k];
        complex double z = in->z;
        complex double c = in->c;
        int index = in->index;

        bool active = true;
        for (int i = 1; i <= chunk; ++i) {

          z = 2 * z_d[i-1] * z + z * z + c;
          complex double rz = z_d[i] + z;
          double rz2 = cnorm(rz);
  
          if (rz2 < z_size[i]) {
            // glitched
            int my_glitch_count;
            #pragma omp atomic capture
            my_glitch_count = glitch_count++;
            struct pixel *out = &glitched[my_glitch_count];
            out->c = c;
            out->z = rz;
            out->index = index;
            out->iters = iters + i;
            active = false;
            break;
          } else if (rz2 > escape_radius_2) {
            // escaped
            img->output[4 * index + 0] = iters + i;
            img->output[4 * index + 1] = 1 - log2(log(rz2) / log_escape_radius_2);
            img->output[4 * index + 2] = carg(rz) / twopi;
            img->output[4 * index + 3] = start_id;
            active = false;
            break;
          }
        }
        if (active) {
          // still iterating
          int my_active_count;
          #pragma omp atomic capture
          my_active_count = active_count++;
          struct pixel *out = &ref->px[1][my_active_count];
          out->c = c;
          out->z = z;
          out->index = index;
        }
      }
  
      if (glitch_count) {
        qsort(glitched, glitch_count, sizeof(*glitched), cmp_pixel_by_iters_asc);
        int end = 0;
        for (int start = 0; start < glitch_count; start = end) {

          // compute minimum of span
          uint32_t start_iters = glitched[start].iters;
          double min_z2 = 1.0 / 0.0;
          uint32_t min_ix = start;
          for (end = start; end < glitch_count && glitched[end].iters == start_iters; ++end) {
            double z2 = cnorm(glitched[end].z);
            if (z2 < min_z2) {
              min_z2 = z2;
              min_ix = end;
            }
          }

          // enqueue a new reference
          int count = end - start;
          struct reference *new_ref = calloc(1, sizeof(*new_ref));
          new_ref->parent = ref;
          new_ref->px[0] = malloc(count * sizeof(*ref->px[0]));
          memcpy(new_ref->px[0], glitched + start, count * sizeof(*ref->px[0]));
          new_ref->count = count;
          new_ref->index = min_ix - start;
          new_ref->period = start_iters;
          new_ref->iters = start_iters;
          image_enqueue(img, new_ref);

        } // for start
      } // if glitch_count
      if (active_count && active_count < ref->count) {
        free(ref->px[0]);
        ref->px[0] = ref->px[1];
        ref->px[1] = malloc(active_count * sizeof(*ref->px[1]));
        free(glitched);
        glitched = malloc(active_count * sizeof(*glitched));
      } else {
        struct pixel *tmp = ref->px[0];
        ref->px[0] = ref->px[1];
        ref->px[1] = tmp;
      }
      ref->count = active_count;
      if (! active_count) {
        break;
      }
    } // for iters

    free(glitched);
  } // while ref

  free(z_d);
  free(z_size);

  image_log(img, LOG_QUEUE, "         END\n");
  return 0;
}
