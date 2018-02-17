#include "special.hpp"
#include "special_internal.hpp"

float std::ellint_3f(float k, float nu, float phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(nu)) return nu;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_3(k, nu, phi); });
}
