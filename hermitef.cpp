#include "special.hpp"
#include "special_internal.hpp"

float std::hermitef(unsigned n, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::hermite(n, x); });
}
