#pragma once

#include <boost/math/special_functions.hpp>
#include "special.hpp"

template<class F> inline
auto _Boost_call(F f) {
    try {
        return f();
    } catch(boost::math::rounding_error&) {
        throw std::domain_error("FIXME: boost::math::rounding_error");
    } catch(boost::math::evaluation_error&) {
        throw std::domain_error("FIXME: boost::math::evaluation_error");
    }
}
