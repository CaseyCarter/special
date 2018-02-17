#include "special.hpp"
#include "special_internal.hpp"

double std::cyl_neumann(double nu, double x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_neumann(nu, x); });
}
