#include "special.hpp"
#include "special_internal.hpp"

double std::ellint_3(double k, double nu, double phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(nu)) return nu;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_3(k, nu, phi); });
}
