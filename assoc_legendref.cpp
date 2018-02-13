#include "special.hpp"
#include "special_internal.hpp"

float std::assoc_legendref(unsigned l, unsigned m, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, m, x); });
}
