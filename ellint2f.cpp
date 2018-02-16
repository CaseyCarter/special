#include "special.hpp"
#include "special_internal.hpp"

float std::ellint_2f(float k, float phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_2(k, phi); });
}
