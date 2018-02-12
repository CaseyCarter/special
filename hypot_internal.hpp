#pragma once

#include <cmath>
#include <limits>
#include <utility>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/policies/error_handling.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/config/no_tr1/cmath.hpp>
#include "special_internal.hpp"

template<class T>
T _Hypot3(T dx, T dy, T dz) _NOEXCEPT
{
    constexpr T eps = boost::math::tools::epsilon<T>();
    constexpr T inf = std::numeric_limits<T>::infinity();

    dx = std::abs(dx);
    dy = std::abs(dy);
    dz = std::abs(dz);

    if (dx == inf || dy == inf || dz == inf) return inf;

    if (dy > dx) std::swap(dx, dy);
    if (dz > dx) std::swap(dx, dz);

    auto fracsq = [](T const numerator, T const denominator) {
        T const result = numerator / denominator;
        return result * result;
    };

    if (dx * eps >= dy && dx * eps >= dz) return dx;

    return dx * std::sqrt(1 + fracsq(dy, dx) + fracsq(dz, dx));
}
