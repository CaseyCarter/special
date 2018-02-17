#include "special.hpp"
#include "special_internal.hpp"

double std::riemann_zeta(double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::zeta(x); });
}
