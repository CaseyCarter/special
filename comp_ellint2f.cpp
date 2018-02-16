#include "special.hpp"
#include "special_internal.hpp"

float std::comp_ellint_2f(float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::ellint_2(x); });
}
