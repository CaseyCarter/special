#include "special.hpp"
#include "special_internal.hpp"

float std::sph_legendref(unsigned l, unsigned m, float theta) {
    if (std::isnan(theta)) return theta;
    return _Boost_call([=]{ return boost::math::spherical_harmonic_r(l, m, theta, 0.0f); });
}
