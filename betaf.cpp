#include "special.hpp"
#include "special_internal.hpp"

float std::betaf(float x, float y) {
    return _Boost_call([=]{ return boost::math::beta(x, y); });
}
