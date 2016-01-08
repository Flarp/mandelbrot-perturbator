// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015,2016 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

template <typename R>
struct pixel {
  std::complex<R> c;
  std::complex<R> z;
  std::complex<R> dz;
  uint32_t index;
  uint32_t iters;
};

template <typename R>
static int cmp_pixel_by_iters_asc(const void *a, const void *b) {
  const struct pixel<R> *x = (const struct pixel<R> *) a;
  const struct pixel<R> *y = (const struct pixel<R> *) b;
  if (x->iters < y->iters) { return -1; }
  if (x->iters > y->iters) { return  1; }
/*
  R x2 = std::norm(x->z);
  R y2 = std::norm(y->z);
  if (x2 < y2) { return -1; }
  if (x2 > y2) { return  1; }
*/
  return 0;
}

template <typename R>
struct reference {
  struct reference<R> *next;

  int has_parent;
  mpc_t parent_c;
  mpc_t parent_z;

  int queue_id;
  int start_id;

  int period;
  int iters;
  mpc_t c;
  mpc_t z;
  std::complex<R> z_d_old;
  std::complex<R> z_d;
  R z_d2eM6;

  struct pixel<R> *px[2];
  int count;

  int index;
};

template <typename R>
static void reference_release(struct reference<R> *ref) {
  if (ref) {
    if (ref->px[0]) {
      free(ref->px[0]);
      ref->px[0] = 0;
    }
    if (ref->px[1]) {
      free(ref->px[1]);
      ref->px[1] = 0;
    }
    if (ref->has_parent) {
      mpc_clear(ref->parent_c);
      mpc_clear(ref->parent_z);
      ref->has_parent = 0;
    }
    mpc_clear(ref->c);
    mpc_clear(ref->z);
    free(ref);
  }
}


template <typename R>
static struct reference<R> *image_dequeue(struct perturbator *img, struct reference<R> *ref) {
  pthread_mutex_lock(&img->mutex);
//  assert(img->ft == FT);
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
    ref = (struct reference<R> *) img->refs;
    img->refs = ref->next;
    ref->start_id = (img->start_id += 1);
    img->active_workers += 1;
    image_log(img, LOG_QUEUE, "%8d START  %8d%8d%8d%16d\n", ref->start_id, ref->queue_id, ref->iters, ref->period, ref->count);
  } else {
    ref = 0;
  }
  pthread_mutex_unlock(&img->mutex);
  return ref;
}

template <typename R>
static void image_enqueue(struct perturbator *img, struct reference<R> *ref) {
  pthread_mutex_lock(&img->mutex);
//  assert(img->ft == FT);
  ref->queue_id = (img->queue_id += 1);
  image_log(img, LOG_QUEUE, "         QUEUE  %8d%8d%8d%16d\n", ref->queue_id, ref->iters, ref->period, ref->count);
  // initialize
  if (ref->has_parent) {
    mpc_init2(ref->z, img->precision);
    mpc_set_ui_ui(ref->z, 0, 0, MPC_RNDNN);
    mpc_init2(ref->c, img->precision);
    mpc_set_ui_ui(ref->c, 0, 0, MPC_RNDNN);
  }
  // link into reference queue (sorted by count descending)
  struct reference<R> *before = 0;
  struct reference<R> *after = (struct reference<R> *) img->refs;
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

template <typename R>
static void *image_worker(void *arg);

template <typename R>
void perturbator_start_internal(struct perturbator *img) {
//  assert(img->ft == FT);
  struct reference<R> *ref = (struct reference<R> *) calloc(1, sizeof(*ref));
  mpc_init2(ref->c, img->precision);
  mpc_init2(ref->z, img->precision);
  mpc_set(ref->c, img->center, MPC_RNDNN);
  mpc_set_ui_ui(ref->z, 0, 0, MPC_RNDNN);
  ref->iters = 0;
  ref->period = 0;
  ref->count = img->width * img->height;
  ref->px[0] = (struct pixel<R> *) calloc(1, ref->count * sizeof(*ref->px[0]));
  R vdiameter = R(2.0) * mpfr_get(img->radius, MPFR_RNDN, R(0));
  R hdiameter = img->width * vdiameter / img->height;
  #pragma omp parallel for
  for (int j = 0; j < img->height; ++j) {
    R y = R((j + 0.5) / img->height - 0.5) * vdiameter;
    for (int i = 0; i < img->width; ++i) {
      R x = R((i + 0.5) / img->width - 0.5) * hdiameter;
      std::complex<R> c(x, -y);
      int index = j * img->width + i;
      ref->px[0][index].c = c;
      ref->px[0][index].z = 0;
      ref->px[0][index].dz = 0;
      ref->px[0][index].index = index;
    }
  }
  img->threads = (pthread_t *) calloc(1, img->workers * sizeof(*img->threads));
  img->running = true;
  image_enqueue(img, ref);
  pthread_mutex_lock(&img->mutex);
  img->active_workers = img->workers;
  for (int i = 0; i < img->workers; ++i) {
    pthread_create(&img->threads[i], 0, image_worker<R>, img);
  }
  pthread_mutex_unlock(&img->mutex);
}

template <typename R>
void release_refs(struct perturbator *img) {
//  assert(img->ft == FT);
  for (struct reference<R> *ref = (struct reference<R> *) img->refs; ref; ) {
    struct reference<R> *next = ref->next;
    reference_release(ref);
    ref = next;
  }
  img->refs = 0;
}

template <typename R>
static void *image_worker(void *arg) {
  struct perturbator *img = (struct perturbator *) arg;
  image_log(img, LOG_QUEUE, "         ENTER\n");

  pthread_mutex_lock(&img->mutex);

  int maxiters = img->maxiters;
  R escape_radius_2 = img->escape_radius_2;
  R log_escape_radius_2 = img->log_escape_radius_2;
  int newton_steps_root = img->newton_steps_root;
//  int newton_steps_child = img->newton_steps_child;
  int chunk = img->chunk;
  R glitch_threshold = img->glitch_threshold;
  int precision = img->precision;
  mpc_t center;
  mpc_init2(center, mpc_get_prec(img->center));
  mpc_set(center, img->center, MPC_RNDNN);
  mpfr_t radius;
  mpfr_init2(radius, 53);
  mpfr_set(radius, img->radius, MPFR_RNDN);
  R pixel_spacing = mpfr_get(radius, MPFR_RNDN, R(0)) * R(2.0 / img->height);
  mpc_t last_reference;
  mpc_init2(last_reference, mpc_get_prec(img->last_reference));
  mpc_set(last_reference, img->last_reference, MPC_RNDNN);
  int last_period = img->last_period;
  float *output = img->output;

  pthread_mutex_unlock(&img->mutex);

  mpc_t *z_hi = (mpc_t *) calloc(1, (1 + chunk) * sizeof(*z_hi));
  for (int k = 0; k <= chunk; ++k) {
    mpc_init2(z_hi[k], precision);
  }
  std::complex<R> *z_d = (std::complex<R> *) calloc(1, (1 + chunk) * sizeof(*z_d));
  R *z_size = (R *) calloc(1, (1 + chunk) * sizeof(*z_size));

  struct reference<R> *ref = 0;
  while ( (ref = image_dequeue(img, ref)) ) {

    // find reference atom
    if (ref->has_parent) {

      // using .z .dz to do one Newton step instead of using nucleus...
      // nucleus = c - z / dz
      mpc_add(ref->c, ref->parent_c, ref->px[0][ref->index].c, MPC_RNDNN);
      mpc_add(ref->z, ref->parent_z, ref->px[0][ref->index].z, MPC_RNDNN);
      mpc_div(ref->z, ref->z,        ref->px[0][ref->index].dz, MPC_RNDNN);
      if (mpfr_number_p(mpc_realref(ref->z)) && mpfr_number_p(mpc_imagref(ref->z))) {
        mpc_sub(ref->c, ref->c, ref->z, MPC_RNDNN);
      }
      // compute rebase offset
      mpc_sub(ref->z, ref->parent_c, ref->c, MPC_RNDNN);
      std::complex<R> dc = -mpc_get(ref->z, MPC_RNDNN, R(0));
      mpc_set_ui_ui(ref->z, 0, 0, MPC_RNDNN);
      // rebase pixels to new reference
      int count = ref->count;
      for (int k = 0; k < count; ++k) {
        ref->px[0][k].c += dc;
      }

    } else {

      int exponent = mpfr_get_exp(radius);
      // find an appropriate initial reference
      bool ok = false;
      bool reused = false;
      int period = 0;
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
        mpfr_sqr(radius2, radius, MPFR_RNDN);
        mpfr_mul_d(radius2, radius2, 65536, MPFR_RNDN);
        mpc_sub(delta, last_reference, center, MPC_RNDNN);
        mpc_norm(delta2, delta, MPFR_RNDN);
        if (mpfr_less_p(delta2, radius2)) {
          mpc_set(nucleus, last_reference, MPC_RNDNN);
          period = last_period;
          if (! mpfr_zero_p(delta2)) {
            exponent = max(exponent, int(mpfr_get_exp(delta2) / 2));
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
        period = m_r_box_period_do(center, radius, maxiters);
        if (period) {
          m_r_nucleus(nucleus, center, period, newton_steps_root);
          if (m_cardioid == m_r_shape(nucleus, period)) {
            ok = true;
            image_log(img, LOG_CACHE, "         BOXED                  %8d\n", period);
          }
        }
      }
      if (0 && ! ok) {
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
        mpfr_sqr(radius2, radius, MPFR_RNDN);
        mpfr_mul_d(radius2, radius2, 65536, MPFR_RNDN);

        for (period = 1; period < maxiters; ++period) {
          if (! image_running(img)) {
            break;
          }
          mpc_sqr(z, z, MPC_RNDNN);
          mpc_add(z, z, center, MPC_RNDNN);
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
        // reuse fallback
        mpc_t delta;
        mpc_init2(delta, 53);
        mpfr_t delta2, radius2;
        mpfr_init2(delta2, 53);
        mpfr_init2(radius2, 53);
        mpfr_sqr(radius2, radius, MPFR_RNDN);
        mpc_sub(delta, last_reference, center, MPC_RNDNN);
        mpc_norm(delta2, delta, MPFR_RNDN);
        mpfr_add(delta2, delta2, radius2, MPFR_RNDN);
        mpc_set(nucleus, last_reference, MPC_RNDNN);
        period = last_period;
        exponent = max(exponent, int(mpfr_get_exp(delta2) / 2));
        ok = true;
        reused = true;
        mpc_clear(delta);
        mpfr_clear(delta2);
        mpfr_clear(radius2);
        image_log(img, LOG_CACHE, "         REUSE FALLBACK\n");
      }

      if (ok) {
        pthread_mutex_lock(&img->mutex);
        mpc_set_prec(img->last_reference, img->precision);
        mpc_set(img->last_reference, nucleus, MPC_RNDNN);
        img->last_period = period;
        pthread_mutex_unlock(&img->mutex);
      }
      mpc_set(ref->c, nucleus, MPC_RNDNN);
      mpc_sub(nucleus, ref->c, center, MPC_RNDNN);
      std::complex<R> dc = -mpc_get(nucleus, MPC_RNDNN, R(0));
      mpc_clear(nucleus);
      mpc_set_ui_ui(ref->z, 0, 0, MPC_RNDNN);
      ref->z_d_old = 0;
      ref->z_d = 0;
      struct series_node *snode = image_cached_approx(img, reused, ref->c, exponent, ref->z, &ref->iters);
      ref->z_d = mpc_get(ref->z, MPC_RNDNN, R(0));

      // rebase pixels to new reference with approximation
      int count = ref->count;
      switch (snode->ft) {
        case ft_float: {
          struct z2c_approx<float> *approx = snode->u.approxf;
          #pragma omp parallel for
          for (int k = 0; k < count; ++k) {
            ref->px[0][k].c += dc;
            std::complex<float> z, dz, c(float(to_ld(std::real(ref->px[0][k].c))), float(to_ld(std::imag(ref->px[0][k].c))));
            z2c_approx_do(approx, c, &z, &dz);
            ref->px[0][k].z = z;
            ref->px[0][k].dz = dz;
          }
          break;
        }
        case ft_double: {
          struct z2c_approx<double> *approx = snode->u.approx;
          #pragma omp parallel for
          for (int k = 0; k < count; ++k) {
            ref->px[0][k].c += dc;
            std::complex<double> z, dz, c(double(to_ld(std::real(ref->px[0][k].c))), double(to_ld(std::imag(ref->px[0][k].c))));
            z2c_approx_do(approx, c, &z, &dz);
            ref->px[0][k].z = z;
            ref->px[0][k].dz = dz;
          }
          break;
        }
        case ft_long_double: {
          struct z2c_approx<long double> *approx = snode->u.approxl;
          #pragma omp parallel for
          for (int k = 0; k < count; ++k) {
            ref->px[0][k].c += dc;
            std::complex<long double> z, dz, c((long double)(to_ld(std::real(ref->px[0][k].c))), (long double)(to_ld(std::imag(ref->px[0][k].c))));
            z2c_approx_do(approx, c, &z, &dz);
            ref->px[0][k].z = z;
            ref->px[0][k].dz = dz;
          }
          break;
        }
        case ft_edouble: {
          struct z2c_approx<edouble> *approx = snode->u.approxe;
          #pragma omp parallel for
          for (int k = 0; k < count; ++k) {
            ref->px[0][k].c += dc;
            std::complex<edouble> z, dz, c(edouble(std::real(ref->px[0][k].c)), edouble(std::imag(ref->px[0][k].c)));
            z2c_approx_do(approx, c, &z, &dz);
            ref->px[0][k].z = to_C(z, R(0));
            ref->px[0][k].dz = to_C(dz, R(0));
          }
          break;
        }
        default: {
          assert(! "valid float type");
        }
      }

      image_log(img, LOG_CACHE, "         APPROX         %8d\n", ref->iters);

    } // if parent

    // prepare for output pixels
    ref->px[1] = (struct pixel<R> *) calloc(1, ref->count * sizeof(*ref->px[1]));
    struct pixel<R> *glitched = (struct pixel<R> *) calloc(1, ref->count * sizeof(*glitched));
//    int start_id = ref->start_id;

    for (int iters = ref->iters; iters < maxiters; iters += chunk) {
      if (! image_running(img)) {
        break;
      }

      // compute ref chunk
      mpc_set(z_hi[0], ref->z, MPC_RNDNN);
      z_d[0] = ref->z_d;
      for (int k = 1; k <= chunk; ++k) {
        mpc_sqr(ref->z, ref->z, MPC_RNDNN);
        mpc_add(ref->z, ref->z, ref->c, MPC_RNDNN);
        mpc_set(z_hi[k], ref->z, MPC_RNDNN);
        z_d[k] = mpc_get(ref->z, MPC_RNDNN, R(0));
      }
      ref->z_d = z_d[chunk];
      for (int k = 0; k <= chunk; ++k) {
        z_size[k] = std::norm(z_d[k]) * glitch_threshold;
      }

      // advance pixels
      int input_count = ref->count;
      int active_count = 0;
      int glitch_count = 0;

      #pragma omp parallel for
      for (int k = 0; k < input_count; ++k) {
        if (image_running(img)) {

          struct pixel<R> *in = &ref->px[0][k];
          std::complex<R> dz = in->dz;
          std::complex<R> z = in->z;
          std::complex<R> c = in->c;
          std::complex<R> rz = z_d[0] + z;
          int index = in->index;

          bool active = true;
          for (int i = 1; i <= chunk; ++i) {

            dz = R(2.0) * rz * dz + R(1.0);
            z = R(2.0) * z_d[i-1] * z + z * z + c;
            rz = z_d[i] + z;
            R rz2 = std::norm(rz);

            if (rz2 < z_size[i]) {
              // glitched
              int my_glitch_count;
              #pragma omp atomic capture
              my_glitch_count = glitch_count++;
              struct pixel<R> *out = &glitched[my_glitch_count];
              out->c = c;
              out->z = rz;
              out->dz = dz;
              out->index = index;
              out->iters = iters + i;
              active = false;
              break;
            } else if (rz2 > escape_radius_2) {
              // escaped
              output[4 * index + 0] = iters + i;
              output[4 * index + 1] = 1 - log2(log(to_ld(rz2)) / to_ld(log_escape_radius_2)); // smooth iters
              output[4 * index + 2] = std::arg(std::complex<long double>(to_ld(std::real(rz)), to_ld(std::imag(rz)))) / twopi;
              output[4 * index + 3] = sqrt(to_ld(rz2)) * log(to_ld(rz2)) / sqrt(to_ld(std::norm(dz * pixel_spacing))); // de
              active = false;
              break;
            }

          } // for i
          if (active) {
            // still iterating
            int my_active_count;
            #pragma omp atomic capture
            my_active_count = active_count++;
            struct pixel<R> *out = &ref->px[1][my_active_count];
            out->c = c;
            out->z = z;
            out->dz = dz;
            out->index = index;
          }

        } // if running
      } // for k

      if (! image_running(img)) { break; }

      if (glitch_count) {
        qsort(glitched, glitch_count, sizeof(*glitched), cmp_pixel_by_iters_asc<R>);
        int end = 0;
        for (int start = 0; start < glitch_count; start = end) {

          // compute minimum of span
          uint32_t start_iters = glitched[start].iters;
          R min_z2 = 1.0 / 0.0;
          uint32_t min_ix = start;
          for (end = start; end < glitch_count && glitched[end].iters == start_iters; ++end) {
            R z2 = std::norm(glitched[end].z);
            if (z2 < min_z2) {
              min_z2 = z2;
              min_ix = end;
            }
          }

          // enqueue a new reference
          int count = end - start;
          struct reference<R> *new_ref = (struct reference<R> *) calloc(1, sizeof(*new_ref));
          // copy parent
          new_ref->has_parent = 1;
          mpc_init2(new_ref->parent_c, img->precision);
          mpc_set(new_ref->parent_c, ref->c, MPC_RNDNN);
          mpc_init2(new_ref->parent_z, img->precision);
          mpc_set(new_ref->parent_z, z_hi[start_iters - iters], MPC_RNDNN);
          new_ref->period = start_iters;
          new_ref->iters = start_iters;
          // copy pixels
          new_ref->px[0] = (struct pixel<R> *) calloc(1, count * sizeof(*ref->px[0]));
          memcpy(new_ref->px[0], glitched + start, count * sizeof(*ref->px[0]));
          new_ref->count = count;
          new_ref->index = min_ix - start;
          image_enqueue(img, new_ref);

        } // for start
      } // if glitch_count
      if (active_count && active_count < ref->count) {
        free(ref->px[0]);
        ref->px[0] = ref->px[1];
        ref->px[1] = (struct pixel<R> *) calloc(1, active_count * sizeof(*ref->px[1]));
        free(glitched);
        glitched = (struct pixel<R> *) calloc(1, active_count * sizeof(*glitched));
      } else {
        struct pixel<R> *tmp = ref->px[0];
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

  for (int k = 0; k <= chunk; ++k) {
    mpc_clear(z_hi[k]);
  }
  free(z_hi);
  free(z_d);
  free(z_size);

  image_log(img, LOG_QUEUE, "         END\n");
  return 0;
}
