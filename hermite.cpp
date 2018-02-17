#include "special.hpp"
#include "special_internal.hpp"

double std::hermite(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::hermite(n, x); });
}
