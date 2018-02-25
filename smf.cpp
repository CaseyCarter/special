#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/math/special_functions.hpp>
#include "special.hpp"

namespace {
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
}

double std::assoc_laguerre(unsigned n, unsigned m, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::laguerre(n, m, x); });
}

float std::assoc_laguerref(unsigned n, unsigned m, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::laguerre(n, m, x); });
}

double std::assoc_legendre(unsigned l, unsigned m, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, m, x); });
}

float std::assoc_legendref(unsigned l, unsigned m, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, m, x); });
}

double std::beta(double x, double y) {
    return _Boost_call([=]{ return boost::math::beta(x, y); });
}

float std::betaf(float x, float y) {
    return _Boost_call([=]{ return boost::math::beta(x, y); });
}

double std::comp_ellint_1(double k) {
    return _Boost_call([=]{ return boost::math::ellint_1(k); });
}

float std::comp_ellint_1f(float k) {
    return _Boost_call([=]{ return boost::math::ellint_1(k); });
}

double std::comp_ellint_2(double k) {
    if (std::isnan(k)) return k;
    return _Boost_call([=]{ return boost::math::ellint_2(k); });
}

float std::comp_ellint_2f(float k) {
    if (std::isnan(k)) return k;
    return _Boost_call([=]{ return boost::math::ellint_2(k); });
}

double std::comp_ellint_3(double k, double nu) {
    if (std::isnan(k)) return k;
    if (std::isnan(nu)) return nu;
    return _Boost_call([=]{ return boost::math::ellint_3(k, nu); });
}

float std::comp_ellint_3f(float k, float nu) {
    if (std::isnan(k)) return k;
    if (std::isnan(nu)) return nu;
    return _Boost_call([=]{ return boost::math::ellint_3(k, nu); });
}

double std::cyl_bessel_i(double nu, double x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_i(nu, x); });
}

float std::cyl_bessel_if(float nu, float x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_i(nu, x); });
}

double std::cyl_bessel_j(double nu, double x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_j(nu, x); });
}

float std::cyl_bessel_jf(float nu, float x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_j(nu, x); });
}

double std::cyl_bessel_k(double nu, double x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_k(nu, x); });
}

float std::cyl_bessel_kf(float nu, float x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_bessel_k(nu, x); });
}

double std::cyl_neumann(double nu, double x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_neumann(nu, x); });
}

float std::cyl_neumannf(float nu, float x) {
    if (std::isnan(nu)) return nu;
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::cyl_neumann(nu, x); });
}

double std::ellint_1(double k, double phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_1(k, phi); });
}

float std::ellint_1f(float k, float phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_1(k, phi); });
}

double std::ellint_2(double k, double phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_2(k, phi); });
}

float std::ellint_2f(float k, float phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_2(k, phi); });
}

double std::ellint_3(double k, double nu, double phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(nu)) return nu;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_3(k, nu, phi); });
}

float std::ellint_3f(float k, float nu, float phi) {
    if (std::isnan(k)) return k;
    if (std::isnan(nu)) return nu;
    if (std::isnan(phi)) return phi;
    return _Boost_call([=]{ return boost::math::ellint_3(k, nu, phi); });
}

double std::expint(double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::expint(x); });
}

float std::expintf(float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::expint(x); });
}

double std::hermite(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::hermite(n, x); });
}

float std::hermitef(unsigned n, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::hermite(n, x); });
}

double std::laguerre(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::laguerre(n, x); });
}

float std::laguerref(unsigned n, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::laguerre(n, x); });
}

double std::legendre(unsigned l, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, x); });
}

float std::legendref(unsigned l, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::legendre_p(l, x); });
}

double std::riemann_zeta(double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::zeta(x); });
}

float std::riemann_zetaf(float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::zeta(x); });
}

double std::sph_bessel(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::sph_bessel(n, x); });
}

float std::sph_besself(unsigned n, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::sph_bessel(n, x); });
}

double std::sph_legendre(unsigned l, unsigned m, double theta) {
    if (std::isnan(theta)) return theta;
    return _Boost_call([=]{ return boost::math::spherical_harmonic_r(l, m, theta, 0.0); });
}

float std::sph_legendref(unsigned l, unsigned m, float theta) {
    if (std::isnan(theta)) return theta;
    return _Boost_call([=]{ return boost::math::spherical_harmonic_r(l, m, theta, 0.0f); });
}

double std::sph_neumann(unsigned n, double x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::sph_neumann(n, x); });
}

float std::sph_neumannf(unsigned n, float x) {
    if (std::isnan(x)) return x;
    return _Boost_call([=]{ return boost::math::sph_neumann(n, x); });
}

namespace {
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
}

double std::hypot(double dx, double dy, double dz) _NOEXCEPT
{
    return _Hypot3<double>(dx, dy, dz);
}

float std::hypot(float dx, float dy, float dz) _NOEXCEPT
{
    return _Hypot3<float>(dx, dy, dz);
}
