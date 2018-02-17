#include "special.hpp"
#include "special_internal.hpp"

double std::sph_bessel(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::sph_bessel(n, x); });
}
