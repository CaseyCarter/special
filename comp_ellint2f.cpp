#include "special.hpp"
#include "special_internal.hpp"

float std::comp_ellint_2f(float k) {
    if (std::isnan(k)) return k;
    return _Boost_call([=]{ return boost::math::ellint_2(k); });
}
