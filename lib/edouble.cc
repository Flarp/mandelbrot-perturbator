// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015,2016 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#ifndef EDOUBLE_CC
#define EDOUBLE_CC 1

#include <cmath>
#include <complex>
#include <algorithm>
#include <iostream>

#include <mpfr.h>

// ldexp is ~2x slower than scaling via table of powers of two
//#define EDOUBLE_USE_LDEXP

inline double sign(double x) {
  return (x > 0) - (x < 0);
}

class edouble {
private:
  static const long maxexponent = (1L << 60) - 2;
  static const long minexponent = 2 - (1L << 60);
  static const long supexponent = maxexponent - 2000;
  static const long infexponent = minexponent + 2080;
  double x;
  long e;
public:
  static const double scaling[129];
  inline edouble() : x(0), e(0) { };
  inline edouble(const edouble &that) : x(that.x), e(that.e) { };
  inline edouble(const double &x0, const long &e0) {
    if (x0 == 0 || std::isnan(x0) || std::isinf(x0)) {
      x = x0;
      e = 0;
    } else {
      int tmp(0);
      double x1(std::frexp(x0, &tmp));
      long e1(tmp);
      e1 += e0;
      if (e0 > supexponent || e1 > maxexponent) {
        x = sign(x0) / 0.0;
        e = 0;
      } else if (e0 < infexponent || e1 < minexponent) {
        x = sign(x0) * 0.0;
        e = 0;
      } else {
        x = x1;
        e = e1;
      }
    }
  };
  inline edouble(const long double &x0, const long &e0) {
    if (x0 == 0 || std::isnan(x0) || std::isinf(x0)) {
      x = x0;
      e = 0;
    } else {
      int tmp(0);
      long double x1(frexp(x0, &tmp));
      long e1(tmp);
      e1 += e0;
      if (e0 > supexponent || e1 > maxexponent) {
        x = sign(x0) / 0.0;
        e = 0;
      } else if (e0 < infexponent || e1 < minexponent) {
        x = sign(x0) * 0.0;
        e = 0;
      } else {
        x = x1;
        e = e1;
      }
    }
  };
  inline edouble(const double &x) : edouble(x, 0) { };
  inline edouble(const long double &x) : edouble(x, 0) { };
  inline explicit edouble(const mpfr_t &a) {
    long e(0);
    double x(mpfr_get_d_2exp(&e, a, MPFR_RNDN));
    *this = edouble(x, e);
  }
  inline edouble(const int &x) : edouble(double(x)) { };
  inline ~edouble() { };
  inline long exponent() const { return e; };
  inline void to_mpfr(mpfr_t &m) {
    mpfr_set_d(m, x, MPFR_RNDN);
    mpfr_mul_2si(m, m, e, MPFR_RNDN);
  };
  inline long double to_ld() {
    int tmp(e);
    if (e > long(tmp)) {
      return sign(x) / 0.0;
    }
    if (e < long(tmp)) {
      return sign(x) * 0.0;
    }
    return std::ldexp((long double) x, tmp);
  }
  friend inline bool operator==(const edouble &a, const edouble &b);
  friend inline bool operator!=(const edouble &a, const edouble &b);
  friend inline int compare(const edouble &a, const edouble &b);
  friend inline int sign(const edouble &a);
  friend inline edouble abs(const edouble &a);
  friend inline edouble sqrt(const edouble &a);
  friend inline edouble operator+(const edouble &a, const edouble &b);
  friend inline edouble operator-(const edouble &a);
  friend inline edouble operator-(const edouble &a, const edouble &b);
  friend inline edouble operator*(const edouble &a, const edouble &b);
  friend inline edouble recip(const edouble &a);
  friend inline edouble operator/(const edouble &a, const edouble &b);
  friend inline bool isnan(const edouble &a);
  friend inline bool isinf(const edouble &a);
  friend inline edouble ldexp(const edouble &a, long e);
  friend inline std::ostream& operator<<(std::ostream& o, const edouble& a);
  friend inline std::istream& operator>>(std::istream& i, edouble& a);
  inline edouble &operator+=(const edouble &a) {
    *this = *this + a;
    return *this;
  };
  inline edouble &operator-=(const edouble &a) {
    *this = *this - a;
    return *this;
  };
  inline edouble &operator*=(const edouble &a) {
    *this = *this * a;
    return *this;
  };
  inline edouble &operator/=(const edouble &a) {
    *this = *this / a;
    return *this;
  };
};
// (0:) . take 128 . map (min 1) . iterate (* 2) $ 2 ^^ negate 128
const double edouble::scaling[129] = {0.0,5.877471754111438e-39,1.1754943508222875e-38,2.350988701644575e-38,4.70197740328915e-38,9.4039548065783e-38,1.88079096131566e-37,3.76158192263132e-37,7.52316384526264e-37,1.504632769052528e-36,3.009265538105056e-36,6.018531076210112e-36,1.2037062152420224e-35,2.407412430484045e-35,4.81482486096809e-35,9.62964972193618e-35,1.925929944387236e-34,3.851859888774472e-34,7.703719777548943e-34,1.5407439555097887e-33,3.0814879110195774e-33,6.162975822039155e-33,1.232595164407831e-32,2.465190328815662e-32,4.930380657631324e-32,9.860761315262648e-32,1.9721522630525295e-31,3.944304526105059e-31,7.888609052210118e-31,1.5777218104420236e-30,3.1554436208840472e-30,6.310887241768095e-30,1.262177448353619e-29,2.524354896707238e-29,5.048709793414476e-29,1.0097419586828951e-28,2.0194839173657902e-28,4.0389678347315804e-28,8.077935669463161e-28,1.6155871338926322e-27,3.2311742677852644e-27,6.462348535570529e-27,1.2924697071141057e-26,2.5849394142282115e-26,5.169878828456423e-26,1.0339757656912846e-25,2.0679515313825692e-25,4.1359030627651384e-25,8.271806125530277e-25,1.6543612251060553e-24,3.308722450212111e-24,6.617444900424222e-24,1.3234889800848443e-23,2.6469779601696886e-23,5.293955920339377e-23,1.0587911840678754e-22,2.117582368135751e-22,4.235164736271502e-22,8.470329472543003e-22,1.6940658945086007e-21,3.3881317890172014e-21,6.776263578034403e-21,1.3552527156068805e-20,2.710505431213761e-20,5.421010862427522e-20,1.0842021724855044e-19,2.168404344971009e-19,4.336808689942018e-19,8.673617379884035e-19,1.734723475976807e-18,3.469446951953614e-18,6.938893903907228e-18,1.3877787807814457e-17,2.7755575615628914e-17,5.551115123125783e-17,1.1102230246251565e-16,2.220446049250313e-16,4.440892098500626e-16,8.881784197001252e-16,1.7763568394002505e-15,3.552713678800501e-15,7.105427357601002e-15,1.4210854715202004e-14,2.842170943040401e-14,5.684341886080802e-14,1.1368683772161603e-13,2.2737367544323206e-13,4.547473508864641e-13,9.094947017729282e-13,1.8189894035458565e-12,3.637978807091713e-12,7.275957614183426e-12,1.4551915228366852e-11,2.9103830456733704e-11,5.820766091346741e-11,1.1641532182693481e-10,2.3283064365386963e-10,4.656612873077393e-10,9.313225746154785e-10,1.862645149230957e-9,3.725290298461914e-9,7.450580596923828e-9,1.4901161193847656e-8,2.9802322387695313e-8,5.960464477539063e-8,1.1920928955078125e-7,2.384185791015625e-7,4.76837158203125e-7,9.5367431640625e-7,1.9073486328125e-6,3.814697265625e-6,7.62939453125e-6,1.52587890625e-5,3.0517578125e-5,6.103515625e-5,1.220703125e-4,2.44140625e-4,4.8828125e-4,9.765625e-4,1.953125e-3,3.90625e-3,7.8125e-3,1.5625e-2,3.125e-2,6.25e-2,0.125,0.25,0.5,1.0};

inline bool operator==(const edouble &a, const edouble &b) {
  return a.e == b.e && a.x == b.x;
}

inline bool operator!=(const edouble &a, const edouble &b) {
  return a.e != b.e || a.x != b.x;
}

inline int compare(const double &a, const double &b) {
  return sign(a - b);
}

inline int compare(const edouble &a, const edouble &b) {
  if (a.x == 0.0) {
    return -sign(b.x);
  }
  if (b.x == 0.0) {
    return sign(a.x);
  }
  long e(std::max(a.e, b.e));
  long da(a.e - e);
  long db(b.e - e);
  int ia(da);
  int ib(db);
  if (long(ia) != da) {
    // a -> 0
    return -sign(b.x);
  }
  if (long(ib) != db) {
    // b -> 0
    return sign(a.x);
  }
#ifdef EDOUBLE_USE_LDEXP
  return compare(std::ldexp(a.x, ia), std::ldexp(b.x, ib));
#else
  double sa = edouble::scaling[std::max(da + 128, 0L)];
  double sb = edouble::scaling[std::max(db + 128, 0L)];
  return compare(a.x * sa, b.x * sb);
#endif
}

inline bool operator<(const edouble &a, const edouble &b) {
  return compare(a, b) < 0;
}

inline bool operator<=(const edouble &a, const edouble &b) {
  return compare(a, b) <= 0;
}

inline bool operator>(const edouble &a, const edouble &b) {
  return compare(a, b) > 0;
}

inline bool operator>=(const edouble &a, const edouble &b) {
  return compare(a, b) >= 0;
}

inline int sign(const edouble &a) {
  return sign(a.x);
}

inline edouble abs(const edouble &a) {
  return { std::abs(a.x), a.e };
}

inline edouble sqrt(const edouble &a) {
  return { std::sqrt((a.e & 1) ? 2.0 * a.x : a.x), (a.e & 1) ? (a.e - 1) / 2 : a.e / 2 };
}

inline edouble operator+(const edouble &a, const edouble &b) {
  if (a.x == 0.0) {
    return b;
  }
  if (b.x == 0.0) {
    return a;
  }
  long e(std::max(a.e, b.e));
  long da(a.e - e);
  long db(b.e - e);
  int ia(da);
  int ib(db);
  if (long(ia) != da) {
    // a -> 0
    return b;
  }
  if (long(ib) != db) {
    // b -> 0
    return a;
  }
#ifdef EDOUBLE_USE_LDEXP
  return edouble(std::ldexp(a.x, ia) + std::ldexp(b.x, ib), e);
#else
  double sa = edouble::scaling[std::max(da + 128, 0L)];
  double sb = edouble::scaling[std::max(db + 128, 0L)];
  return edouble(a.x * sa + b.x * sb, e);
#endif
}

inline edouble operator-(const edouble &a) {
  return { -a.x, a.e };
}

inline edouble operator-(const edouble &a, const edouble &b) {
  if (a.x == 0.0) {
    return -b;
  }
  if (b.x == 0.0) {
    return a;
  }
  long e(std::max(a.e, b.e));
  long da(a.e - e);
  long db(b.e - e);
  int ia(da);
  int ib(db);
  if (long(ia) != da) {
    // a -> 0
    return -b;
  }
  if (long(ib) != db) {
    // b -> 0
    return a;
  }
#ifdef EDOUBLE_USE_LDEXP
  return edouble(std::ldexp(a.x, ia) - std::ldexp(b.x, ib), e);
#else
  double sa = edouble::scaling[std::max(da + 128, 0L)];
  double sb = edouble::scaling[std::max(db + 128, 0L)];
  return edouble(a.x * sa - b.x * sb, e);
#endif
}

inline edouble operator*(const edouble &a, const edouble &b) {
  return edouble(a.x * b.x, a.e + b.e);
}

inline double recip(const double &a) {
  return 1.0 / a;
}

inline edouble recip(const edouble &a) {
  return edouble(recip(a.x), -a.e);
}

inline edouble operator/(const edouble &a, const edouble &b) {
  return edouble(a.x / b.x, a.e - b.e);
}

inline bool isnan(const edouble &a) {
  return isnan(a.x);
}

inline bool isinf(const edouble &a) {
  return isinf(a.x);
}

inline edouble ldexp(const edouble &a, long e) {
  return edouble(a.x, a.e + e);
}

inline void to_mpfr(float from, mpfr_t &to) {
  mpfr_set_flt(to, from, MPFR_RNDN);
}

inline void to_mpfr(double from, mpfr_t &to) {
  mpfr_set_d(to, from, MPFR_RNDN);
}

inline void to_mpfr(long double from, mpfr_t &to) {
  mpfr_set_ld(to, from, MPFR_RNDN);
}

inline void to_mpfr(edouble from, mpfr_t &to) {
  from.to_mpfr(to);
}

inline long double to_ld(float x) {
  return x;
}

inline long double to_ld(double x) {
  return x;
}

inline long double to_ld(long double x) {
  return x;
}

inline long double to_ld(edouble x) {
  return x.to_ld();
}

template <typename S, typename T>
inline T to_R(S x, const T &dummy) {
  (void) dummy;
  return T(to_ld(x));
}

inline edouble to_R(edouble x, const edouble &dummy) {
  (void) dummy;
  return x;
}

template <typename S, typename T>
inline std::complex<T> to_C(std::complex<S> x, const T &dummy) {
  return std::complex<T>(to_R(std::real(x), dummy), to_R(std::imag(x), dummy));
}

inline std::complex<edouble> to_C(std::complex<edouble> x, const edouble &dummy) {
  (void) dummy;
  return x;
}

inline std::ostream& operator<<(std::ostream& o, const edouble& a) {
  return o << a.x << " " << a.e;
}

inline std::istream& operator>>(std::istream& i, edouble& a) {
  double x;
  long e;
  i >> x >> e;
  a = edouble(x, e);
  return i;
}

long exponent(float z)
{
  int e;
  frexp(z, &e);
  return e;
}

long exponent(double z)
{
  int e;
  frexp(z, &e);
  return e;
}

long exponent(long double z)
{
  int e;
  frexpl(z, &e);
  return e;
}

long exponent(edouble z)
{
  return z.exponent();
}

long exponent(const mpfr_t z)
{
  if (mpfr_regular_p(z))
    return mpfr_get_exp(z);
  return 0;
}

bool isfinite(edouble z)
{
  return !(isinf(z) || isnan(z));
}

#ifdef EDOUBLE_STANDALONE

#include <cstdlib>
#include <cstdio>

edouble erandom()
{
  int r = std::rand();
  return edouble((1.0 + (r / (double) RAND_MAX)) * ((r & 1) - 0.5), (r & 511) - 256);
}

int main() {
#if 0
  mpfr_t m;
  mpfr_init2(m, 53);
  edouble x(2);
  do {
    std::cout << x << std::endl;
    x.to_mpfr(m);
    mpfr_printf("%Re\n", m);
    x *= x;
  } while (! isinf(x));
  mpfr_clear(m);
#endif
#define N (1<<14)
  edouble r[N];
  for (int i = 0; i < N; ++i)
    r[i] = erandom();
  printf("+\n"); fflush(stdout);
  edouble sum(0);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
      sum += r[i] + r[j];
  std::cerr << sum << std::endl;
  printf("DONE\n"); fflush(stdout);
  return 0;
}

#endif

#endif
