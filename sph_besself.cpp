#include "special.hpp"
#include "special_internal.hpp"

float std::sph_besself(unsigned n, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::sph_bessel(n, x); });
}