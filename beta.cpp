#include <special.hpp>
#include "special_internal.hpp"

double std::beta(double x, double y) {
    return _Boost_call([=]{ return boost::math::beta(x, y); });
}
