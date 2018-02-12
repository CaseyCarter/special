#include "hypot_internal.hpp"

double std::hypot(double dx, double dy, double dz) _NOEXCEPT
{
    return _Hypot3<double>(dx, dy, dz);
}
