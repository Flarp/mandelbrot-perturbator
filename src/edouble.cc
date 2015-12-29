#ifndef EDOUBLE_CC
#define EDOUBLE_CC 1

#include <cmath>
#include <algorithm>
#include <iostream>

#include <mpfr.h>

double sign(double x) {
  return (x > 0) - (x < 0);
}

class edouble {
private:
  static const int maxexponent = (1 << 30) - 2;
  static const int minexponent = 2 - (1 << 30);
  static const int supexponent = maxexponent - 2000;
  static const int infexponent = minexponent + 2080;
  double x;
  int e;
public:
  edouble() : x(0), e(0) { };
  edouble(const edouble &that) : x(that.x), e(that.e) { };
  edouble(const double &x0, const int &e0) {
    if (x0 == 0 || std::isnan(x0) || std::isinf(x0)) {
      x = x0;
      e = 0;
    } else {
      int e1(0);
      double x1(frexp(x0, &e1));
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
  edouble(const double &x) : edouble(x, 0) { };
  edouble(const mpfr_t &a) {
    long e(0);
    double x(mpfr_get_d_2exp(&e, a, MPFR_RNDN));
    *this = edouble(x, e);
  }
  ~edouble() { };
  void to_mpfr(mpfr_t &m) {
    mpfr_set_d(m, x, MPFR_RNDN);
    mpfr_mul_2si(m, m, e, MPFR_RNDN);
  };
  friend bool operator==(const edouble &a, const edouble &b);
  friend bool operator!=(const edouble &a, const edouble &b);
  friend int compare(const edouble &a, const edouble &b);
  friend int sign(const edouble &a);
  friend edouble abs(const edouble &a);
  friend edouble operator+(const edouble &a, const edouble &b);
  friend edouble operator-(const edouble &a);
  friend edouble operator-(const edouble &a, const edouble &b);
  friend edouble operator*(const edouble &a, const edouble &b);
  friend edouble recip(const edouble &a);
  friend edouble operator/(const edouble &a, const edouble &b);
  friend bool isnan(const edouble &a);
  friend bool isinf(const edouble &a);
  friend std::ostream& operator<<(std::ostream& o, const edouble& a);
  friend std::istream& operator>>(std::istream& i, edouble& a);
  edouble &operator+=(const edouble &a) {
    *this = *this + a;
    return *this;
  };
  edouble &operator-=(const edouble &a) {
    *this = *this - a;
    return *this;
  };
  edouble &operator*=(const edouble &a) {
    *this = *this * a;
    return *this;
  };
  edouble &operator/=(const edouble &a) {
    *this = *this / a;
    return *this;
  };
};

bool operator==(const edouble &a, const edouble &b) {
  return a.e == b.e && a.x == b.x;
}

bool operator!=(const edouble &a, const edouble &b) {
  return a.e != b.e || a.x != b.x;
}

int compare(const edouble &a, const edouble &b) {
  int e(std::max(a.e, b.e));
  return compare(ldexp(a.x, a.e - e), ldexp(b.x, b.e - e));
}

bool operator<(const edouble &a, const edouble &b) {
  return compare(a, b) < 0;
}

bool operator<=(const edouble &a, const edouble &b) {
  return compare(a, b) <= 0;
}

bool operator>(const edouble &a, const edouble &b) {
  return compare(a, b) > 0;
}

bool operator>=(const edouble &a, const edouble &b) {
  return compare(a, b) >= 0;
}

int sign(const edouble &a) {
  return sign(a.x);
}

edouble abs(const edouble &a) {
  return { std::abs(a.x), a.e };
}

edouble operator+(const edouble &a, const edouble &b) {
  int e(std::max(a.e, b.e));
  return edouble(ldexp(a.x, a.e - e) + ldexp(b.x, b.e - e), e);
}

edouble operator-(const edouble &a) {
  return { -a.x, a.e };
}

edouble operator-(const edouble &a, const edouble &b) {
  int e(std::max(a.e, b.e));
  return edouble(ldexp(a.x, a.e - e) - ldexp(b.x, b.e - e), e);
}

edouble operator*(const edouble &a, const edouble &b) {
  return edouble(a.x * b.x, a.e + b.e);
}

double recip(const double &a) {
  return 1.0 / a;
}

edouble recip(const edouble &a) {
  return edouble(recip(a.x), -a.e);
}

edouble operator/(const edouble &a, const edouble &b) {
  return edouble(a.x / b.x, a.e - b.e);
}

bool isnan(const edouble &a) {
  return isnan(a.x);
}

bool isinf(const edouble &a) {
  return isinf(a.x);
}

std::ostream& operator<<(std::ostream& o, const edouble& a) {
  return o << a.x << " " << a.e;
}

std::istream& operator>>(std::istream& i, edouble& a) {
  double x;
  int e;
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
