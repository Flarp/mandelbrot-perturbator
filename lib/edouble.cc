// perturbator -- efficient deep zooming for Mandelbrot sets
// Copyright (C) 2015,2016 Claude Heiland-Allen
// License GPL3+ http://www.gnu.org/licenses/gpl.html

#ifndef EDOUBLE_CC
#define EDOUBLE_CC 1

#include <cmath>
#include <algorithm>
#include <iostream>

#include <mpfr.h>

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
  inline long exponent() { return e; } const;
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
  return compare(std::ldexp(a.x, ia), std::ldexp(b.x, ib));
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
  return edouble(std::ldexp(a.x, ia) + std::ldexp(b.x, ib), e);
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
  return edouble(std::ldexp(a.x, ia) - std::ldexp(b.x, ib), e);
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

#ifdef EDOUBLE_STANDALONE

int main() {
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
  return 0;
}

#endif

#endif
