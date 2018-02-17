#include "special.hpp"
#include "special_internal.hpp"

double std::cyl_bessel_j(double nu, double x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_j(nu, x); });
}
