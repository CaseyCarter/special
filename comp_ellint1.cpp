#include "special.hpp"
#include "special_internal.hpp"

double std::comp_ellint_1(double x) {
    return _Boost_call([=]{ return boost::math::ellint_1(x); });
}
