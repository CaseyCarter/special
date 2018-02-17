#include "special.hpp"
#include "special_internal.hpp"

double std::comp_ellint_2(double k) {
    if (std::isnan(k)) return k;
    return _Boost_call([=]{ return boost::math::ellint_2(k); });
}
