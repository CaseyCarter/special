#include "special.hpp"
#include "special_internal.hpp"

float std::cyl_bessel_if(float nu, float x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_i(nu, x); });
}
