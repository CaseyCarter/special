#include <special.hpp>
#include "special_internal.hpp"

float std::legendref(unsigned l, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, x); });
}
