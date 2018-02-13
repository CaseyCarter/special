#pragma once

#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include "special_internal.hpp"

template<class T> inline
T _Hypot3(T dx, T dy, T dz) _NOEXCEPT
{
    static_assert(std::is_floating_point_v<T>);
    dx = std::abs(dx);
    dy = std::abs(dy);
    dz = std::abs(dz);

    constexpr T inf = std::numeric_limits<T>::infinity();
    if (dx == inf || dy == inf || dz == inf) return inf;

    if (dy > dx) std::swap(dx, dy);
    if (dz > dx) std::swap(dx, dz);

    constexpr T eps = boost::math::tools::epsilon<T>();
    if (dx * eps >= dy && dx * eps >= dz) return dx;

    auto fracsq = [](T const numerator, T const denominator) {
        T const result = numerator / denominator;
        return result * result;
    };

    return dx * std::sqrt(1 + fracsq(dy, dx) + fracsq(dz, dx));
}
