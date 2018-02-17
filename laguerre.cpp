#include "special.hpp"
#include "special_internal.hpp"

double std::laguerre(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::laguerre(n, x); });
}
