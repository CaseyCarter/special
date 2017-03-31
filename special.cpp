#include "special.hpp"
#include <boost/math/special_functions.hpp>

_STD_BEGIN

double assoc_laguerre(unsigned n, unsigned m, double x) {
    if (_CSTD isnan(x)) {
        return x;
    } else {
        return boost::math::laguerre(n, m, x);
    }
}

double assoc_legendre(unsigned l, unsigned m, double x) {
    if (_CSTD isnan(x)) {
        return x;
    } else {
        return boost::math::legendre_p(l, m, x);
    }
}

double legendre(unsigned l, double x) {
    if (_CSTD isnan(x)) {
        return x;
    } else {
        return boost::math::legendre_p(l, x);
    }
}

_STD_END
