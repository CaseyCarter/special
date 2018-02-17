#include "special.hpp"
#include "special_internal.hpp"

double std::sph_legendre(unsigned l, unsigned m, double theta) {
    if (std::isnan(theta)) return theta;
    return _Boost_call([=]{ return boost::math::spherical_harmonic_r(l, m, theta, 0.0); });
}
