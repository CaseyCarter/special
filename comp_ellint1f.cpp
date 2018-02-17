#include "special.hpp"
#include "special_internal.hpp"

float std::comp_ellint_1f(float k) {
    return _Boost_call([=]{ return boost::math::ellint_1(k); });
}
