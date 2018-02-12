#include "hypot_internal.hpp"

float std::hypot(float dx, float dy, float dz) _NOEXCEPT
{
    return _Hypot3<float>(dx, dy, dz);
}
