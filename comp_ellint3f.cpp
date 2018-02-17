#include "special.hpp"
#include "special_internal.hpp"

float std::comp_ellint_3f(float k, float nu) {
    if (std::isnan(k)) return k;
    if (std::isnan(nu)) return nu;
    return _Boost_call([=]{ return boost::math::ellint_3(k, nu); });
}
