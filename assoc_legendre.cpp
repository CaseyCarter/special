#include <special.hpp>
#include "special_internal.hpp"

double std::assoc_legendre(unsigned l, unsigned m, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, m, x); });
}
