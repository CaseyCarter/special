#include "special.hpp"
#include "special_internal.hpp"

double std::sph_neumann(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::sph_neumann(n, x); });
}
