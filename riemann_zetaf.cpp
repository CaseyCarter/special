#include "special.hpp"
#include "special_internal.hpp"

float std::riemann_zetaf(float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::zeta(x); });
}
