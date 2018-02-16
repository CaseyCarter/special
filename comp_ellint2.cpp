#include "special.hpp"
#include "special_internal.hpp"

double std::comp_ellint_2(double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::ellint_2(x); });
}
