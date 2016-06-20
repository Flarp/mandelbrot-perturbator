// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015,2016 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#include <complex>
#include <algorithm>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpfr.h>

#include "edouble.cc"

template< typename R >
struct z2c_series_t {
  int order;
  mpfr_t v[11];
  std::complex<R> *a[2];
  std::complex<R> *b[2];
  int n;
};

struct z2c_series {
  float_type tag;
  union {
    z2c_series_t<double> *d;
    z2c_series_t<long double> *l;
    z2c_series_t<edouble> *e;
  } u;
};

template< typename R >
struct z2c_approx_t {
  int order;
  std::complex<R> *a;
  std::complex<R> *b;
};

struct z2c_approx {
  float_type tag;
  union {
    z2c_approx_t<double> *d;
    z2c_approx_t<long double> *l;
    z2c_approx_t<edouble> *e;
  } u;
};

struct z2c_reference {
  mpfr_t v[6];
  int n;
};

struct z2c_reference *z2c_reference_new(const mpfr_t cx, const mpfr_t cy, const mpfr_t zx, const mpfr_t zy, int n) {
  struct z2c_reference *r = (struct z2c_reference *) malloc(sizeof(*r));
  if (! r) { return 0; }
  mpfr_prec_t p = std::max(std::max(mpfr_get_prec(cx), mpfr_get_prec(cy)), std::max(mpfr_get_prec(zx), mpfr_get_prec(zy)));
  for (int i = 0; i < 6; ++i) {
    mpfr_init2(r->v[i], p);
    mpfr_set_si(r->v[i], 0, MPFR_RNDN);
  };
  mpfr_set(r->v[0], cx, MPFR_RNDN);
  mpfr_set(r->v[1], cy, MPFR_RNDN);
  mpfr_set(r->v[2], zx, MPFR_RNDN);
  mpfr_set(r->v[3], zy, MPFR_RNDN);
  r->n = n;
  return r;
}

void z2c_reference_delete(struct z2c_reference *r) {
  for (int i = 0; i < 6; ++i) {
    mpfr_clear(r->v[i]);
  }
  free(r);
}

void z2c_reference_step(struct z2c_reference *r) {
  mpfr_sqr(r->v[4], r->v[2], MPFR_RNDN);
  mpfr_sqr(r->v[5], r->v[3], MPFR_RNDN);
  mpfr_sub(r->v[4], r->v[4], r->v[5], MPFR_RNDN);
  mpfr_mul(r->v[5], r->v[2], r->v[3], MPFR_RNDN);
  mpfr_mul_2si(r->v[5], r->v[5], 1, MPFR_RNDN);
  mpfr_add(r->v[2], r->v[4], r->v[0], MPFR_RNDN);
  mpfr_add(r->v[3], r->v[5], r->v[1], MPFR_RNDN);
  r->n += 1;
}

int z2c_reference_get_n(const struct z2c_reference *r) {
  return r->n;
}

/*
std::complex<float> z2c_reference_get_zf(const struct z2c_reference *r) {
  return std::complex<float>(mpfr_get_flt(r->v[2], MPFR_RNDN), mpfr_get_flt(r->v[3], MPFR_RNDN));
}

std::complex<double> z2c_reference_get_z(const struct z2c_reference *r) {
  return std::complex<double>(mpfr_get_d(r->v[2], MPFR_RNDN), mpfr_get_d(r->v[3], MPFR_RNDN));
}

std::complex<long double> z2c_reference_get_zl(const struct z2c_reference *r) {
  return std::complex<long double>(mpfr_get_ld(r->v[2], MPFR_RNDN), mpfr_get_ld(r->v[3], MPFR_RNDN));
}
*/
void z2c_reference_get_zr(const struct z2c_reference *ref, mpfr_t zx, mpfr_t zy) {
  mpfr_set_prec(zx, mpfr_get_prec(ref->v[2]));
  mpfr_set_prec(zy, mpfr_get_prec(ref->v[3]));
  mpfr_set(zx, ref->v[2], MPFR_RNDN);
  mpfr_set(zy, ref->v[3], MPFR_RNDN);
}

template< typename R >
struct z2c_approx_t<R> *z2c_approx_t_new(const struct z2c_series_t<R> *s) {
  struct z2c_approx_t<R> *a = (struct z2c_approx_t<R> *) malloc(sizeof(*a));
  a->order = s->order;
  a->a = (std::complex<R> *) calloc(1, s->order * sizeof(std::complex<R>));
  a->b = (std::complex<R> *) calloc(1, (s->order + 1) * sizeof(std::complex<R>));
  R dre = mpfr_get(s->v[4], MPFR_RNDN, R(0));
  R dim = mpfr_get(s->v[5], MPFR_RNDN, R(0));
  a->b[0] = std::complex<R>(dre, dim);
  for (int i = 0; i < s->order; ++i) {
    a->a[i] = s->a[0][i];
    a->b[i+1] = s->b[0][i];
  }
  return a;
}

struct z2c_approx *z2c_approx_new(const struct z2c_series *s) {
  struct z2c_approx *a = (struct z2c_approx *) malloc(sizeof(*a));
  switch (s->tag) {
    case ft_double:      a->u.d = z2c_approx_t_new(s->u.d); break;
    case ft_long_double: a->u.l = z2c_approx_t_new(s->u.l); break;
    case ft_edouble:     a->u.e = z2c_approx_t_new(s->u.e); break;
    default: assert(! "float type valid");
  }
  a->tag = s->tag;
  return a;
}

template< typename R >
void z2c_approx_t_delete(struct z2c_approx_t<R> *a) {
  free(a->a);
  free(a->b);
  free(a);
}

void z2c_approx_delete(struct z2c_approx *a) {
  switch (a->tag) {
    case ft_double:      z2c_approx_t_delete(a->u.d); break;
    case ft_long_double: z2c_approx_t_delete(a->u.l); break;
    case ft_edouble:     z2c_approx_t_delete(a->u.e); break;
    default: assert(! "float type valid");
  }
  free(a);
}

template< typename R >
void z2c_approx_t_do(const struct z2c_approx_t<R> *a, std::complex<R> dc, std::complex<R> *dz_out, std::complex<R> *ddz_out) {
  std::complex<R> z(dc);
  std::complex<R> zi(z);
  std::complex<R> s(0);
  std::complex<R> ds(a->b[0]);
  for (int i = 0; i < a->order; ++i) {
    s += a->a[i] * zi;
    ds += a->b[i + 1] * zi;
    zi *= z;
  }
  *dz_out = s;
  *ddz_out = ds;
}

void z2c_approx_do(const struct z2c_approx *a, std::complex<edouble> dc, std::complex<edouble> *dz_out, std::complex<edouble> *ddz_out) {
  switch (a->tag) {
    case ft_double:
    {
      std::complex<double> dc_ = to_C(dc, double(0));
      std::complex<double> dz_out_, ddz_out_;
      z2c_approx_t_do(a->u.d, dc_, &dz_out_, &ddz_out_);
      *dz_out = to_C(dz_out_, edouble(0));
      *ddz_out = to_C(ddz_out_, edouble(0));
      break;
    }
    case ft_long_double:
    {
      std::complex<long double> dc_ = to_C(dc, (long double)(0));
      std::complex<long double> dz_out_, ddz_out_;
      z2c_approx_t_do(a->u.l, dc_, &dz_out_, &ddz_out_);
      *dz_out = to_C(dz_out_, edouble(0));
      *ddz_out = to_C(ddz_out_, edouble(0));
      break;
    }
    case ft_edouble:
    {
      z2c_approx_t_do(a->u.e, dc, dz_out, ddz_out);
      break;
    }
    default:
    {
      assert(! "float type valid");
    }
  }
}

template< typename R >
struct z2c_series_t<R> *z2c_series_t_new(int order, const mpfr_t cx, const mpfr_t cy, const R dummy) {
  (void) dummy;
  struct z2c_series_t<R> *s = (struct z2c_series_t<R> *) malloc(sizeof(*s));
  if (! s) { return 0; }
  s->order = order;
  s->a[0] = (std::complex<R> *) calloc(1, order * sizeof(std::complex<R>));
  s->a[1] = (std::complex<R> *) calloc(1, order * sizeof(std::complex<R>));
  s->b[0] = (std::complex<R> *) calloc(1, order * sizeof(std::complex<R>));
  s->b[1] = (std::complex<R> *) calloc(1, order * sizeof(std::complex<R>));
  mpfr_prec_t p = std::max(mpfr_get_prec(cx), mpfr_get_prec(cy));
  for (int i = 0; i < 11; ++i) {
    mpfr_init2(s->v[i], p);
    mpfr_set_si(s->v[i], 0, MPFR_RNDN);
  };
  mpfr_set(s->v[0], cx, MPFR_RNDN);
  mpfr_set(s->v[1], cy, MPFR_RNDN);
  mpfr_set(s->v[2], cx, MPFR_RNDN);
  mpfr_set(s->v[3], cy, MPFR_RNDN);
  mpfr_set_si(s->v[4], 1, MPFR_RNDN);
  s->a[0][0] = std::complex<R>(1);
  s->n = 1;
  return s;
}

struct z2c_series *z2c_series_new(int order, const mpfr_t cx, const mpfr_t cy, float_type type) {
  struct z2c_series *s = (struct z2c_series *) malloc(sizeof(*s));
  switch (type) {
    case ft_double:      s->u.d = z2c_series_t_new(order, cx, cy, double(0)); break;
    case ft_long_double: s->u.l = z2c_series_t_new(order, cx, cy, (long double)(0)); break;
    case ft_edouble:     s->u.e = z2c_series_t_new(order, cx, cy, edouble(0)); break;
    default: assert(! "float type valid");
  }
  s->tag = type;
  return s;
}

template< typename R >
void z2c_series_t_delete(struct z2c_series_t<R> *s) {
  for (int i = 0; i < 11; ++i) {
    mpfr_clear(s->v[i]);
  }
  free(s->a[0]);
  free(s->a[1]);
  free(s->b[0]);
  free(s->b[1]);
  free(s);
}

void z2c_series_delete(struct z2c_series *s) {
  switch (s->tag) {
    case ft_double:      z2c_series_t_delete(s->u.d); break;
    case ft_long_double: z2c_series_t_delete(s->u.l); break;
    case ft_edouble:     z2c_series_t_delete(s->u.e); break;
    default: assert(! "float type valid");
  }
  free(s);
}

template< typename R >
int z2c_series_t_get_n(const struct z2c_series_t<R> *s) {
  return s->n;
}

int z2c_series_get_n(struct z2c_series *s) {
  switch (s->tag) {
    case ft_double:      return z2c_series_t_get_n(s->u.d);
    case ft_long_double: return z2c_series_t_get_n(s->u.l);
    case ft_edouble:     return z2c_series_t_get_n(s->u.e);
    default: assert(! "float type valid");
  }
  return 0;
}

template< typename R >
bool z2c_series_t_step(struct z2c_series_t<R> *s, mpfr_exp_t exponent1, mpfr_exp_t threshold, int approx_skip) {
  // 6,7 <- z^2 + c
  mpfr_sqr(s->v[6], s->v[2], MPFR_RNDN);
  mpfr_sqr(s->v[7], s->v[3], MPFR_RNDN);
  mpfr_sub(s->v[6], s->v[6], s->v[7], MPFR_RNDN);
  mpfr_mul(s->v[7], s->v[2], s->v[3], MPFR_RNDN);
  mpfr_mul_2si(s->v[7], s->v[7], 1, MPFR_RNDN);
  mpfr_add(s->v[6], s->v[6], s->v[0], MPFR_RNDN);
  mpfr_add(s->v[7], s->v[7], s->v[1], MPFR_RNDN);
  // 8,9 <- 2 z dz + 1
  mpfr_mul(s->v[8], s->v[2], s->v[4], MPFR_RNDN);
  mpfr_mul(s->v[9], s->v[3], s->v[5], MPFR_RNDN);
  mpfr_sub(s->v[8], s->v[8], s->v[9], MPFR_RNDN);
  mpfr_mul_2si(s->v[8], s->v[8], 1, MPFR_RNDN);
  mpfr_add_ui(s->v[8], s->v[8], 1, MPFR_RNDN);
  mpfr_mul(s->v[9], s->v[2], s->v[5], MPFR_RNDN);
  mpfr_mul(s->v[10], s->v[3], s->v[4], MPFR_RNDN);
  mpfr_add(s->v[9], s->v[9], s->v[10], MPFR_RNDN);
  mpfr_mul_2si(s->v[9], s->v[9], 1, MPFR_RNDN);
  // step coefficients
  std::complex<R> z(mpfr_get(s->v[2], MPFR_RNDN, R(0)), mpfr_get(s->v[3], MPFR_RNDN, R(0)));
  std::complex<R> dz(mpfr_get(s->v[4], MPFR_RNDN, R(0)), mpfr_get(s->v[5], MPFR_RNDN, R(0)));
  R two(2);
  #pragma omp parallel for
  for (int i = 0; i < s->order/2; ++i) {
    // load balancing
    int ns[2] = { i + 1, s->order - i };
    for (int j = 0; j < 2; ++j) {
      int n = ns[j];
      // special case
      if (n == 1) {
        s->a[1][n-1] = two * z * s->a[0][n-1] + R(1);
      } else {
        // a(n) <- 2 z a(n) + sum a(m) a(n-m)
        std::complex<R> sum(0);
        for (int m = 1; m < (n+1)/2; ++m) {
          sum += s->a[0][m - 1] * s->a[0][n-m - 1];
        }
        sum *= two;
        if (0 == (n & 1)) {
          sum += s->a[0][n/2 - 1] * s->a[0][n/2 - 1];
        }
        s->a[1][n-1] = two * z * s->a[0][n-1] + sum;
      }
      // b(n) <- 2 (z b(n) + dz a(n) + sum a(m) b(n-m) )
      std::complex<R> sum(0);
      for (int m = 1; m < n; ++m) {
        sum += s->a[0][m - 1] * s->b[0][n-m - 1];
      }
      s->b[1][n-1] = two * (z * s->b[0][n-1] + dz * s->a[0][n-1] + sum);
    }
  }
  // check valid
  bool valid = true, dvalid = true;
  mpfr_exp_t e0;
  mpfr_exp_t e1;
  mpfr_exp_t de0;
  mpfr_exp_t de1;
  for (int i = 0; i < s->order; ++i) {
    e1 = LONG_MIN;
    de1 = LONG_MIN;
    if (real(s->a[1][i]) != 0) { e1 = std::max(e1, exponent(real(s->a[1][i]))); }
    if (imag(s->a[1][i]) != 0) { e1 = std::max(e1, exponent(imag(s->a[1][i]))); }
    if (real(s->b[1][i]) != 0) { de1 = std::max(de1, exponent(real(s->b[1][i]))); }
    if (imag(s->b[1][i]) != 0) { de1 = std::max(de1, exponent(imag(s->b[1][i]))); }
    if (i > 0) {
      valid = e0 - exponent1 >= e1 + threshold;
      dvalid = de0 - exponent1 >= de1 + threshold;
      e0 = std::max(e0 - exponent1, e1);
      de0 = std::max(de0 - exponent1, de1);
    } else {
      e0 = e1;
      de0 = de1;
    }
  }
  if (approx_skip) {
    valid = dvalid = s->n < approx_skip;
  }
  if ((! valid) || (! dvalid)) { return false; }
  // step
  s->n += 1;
  for (int i = 0; i < 4; ++i) {
    mpfr_set(s->v[i+2], s->v[i+6], MPFR_RNDN);
  }
  for (int i = 0; i < s->order; ++i) {
    s->a[0][i] = s->a[1][i];
  }
  for (int i = 0; i < s->order; ++i) {
    s->b[0][i] = s->b[1][i];
  }
  // ok
  return true;
}

bool z2c_series_step(struct z2c_series *s, mpfr_exp_t exponent1, mpfr_exp_t threshold, int approx_skip) {
  switch (s->tag) {
    case ft_double:      return z2c_series_t_step(s->u.d, exponent1, threshold, approx_skip);
    case ft_long_double: return z2c_series_t_step(s->u.l, exponent1, threshold, approx_skip);
    case ft_edouble:     return z2c_series_t_step(s->u.e, exponent1, threshold, approx_skip);
    default: assert(! "float type valid");
  }
  return 0;
}

template< typename R >
struct z2c_reference *z2c_series_t_reference_new(const struct z2c_series_t<R> *s) {
  return z2c_reference_new(s->v[0], s->v[1], s->v[2], s->v[3], s->n);
}

struct z2c_reference *z2c_series_reference_new(const struct z2c_series *s) {
  switch (s->tag) {
    case ft_double:      return z2c_series_t_reference_new(s->u.d);
    case ft_long_double: return z2c_series_t_reference_new(s->u.l);
    case ft_edouble:     return z2c_series_t_reference_new(s->u.e);
    default: assert(! "float type valid");
  }
  return 0;
}
