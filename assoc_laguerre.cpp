#include "special.hpp"
#include "special_internal.hpp"

double std::assoc_laguerre(unsigned n, unsigned m, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::laguerre(n, m, x); });
}
