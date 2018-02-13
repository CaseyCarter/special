#include "special.hpp"
#include "special_internal.hpp"

float std::assoc_laguerref(unsigned n, unsigned m, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::laguerre(n, m, x); });
}
