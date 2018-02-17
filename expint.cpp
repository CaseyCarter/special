#include "special.hpp"
#include "special_internal.hpp"

double std::expint(double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::expint(x); });
}
