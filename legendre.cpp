#include <special.hpp>
#include "special_internal.hpp"

double std::legendre(unsigned l, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, x); });
}
