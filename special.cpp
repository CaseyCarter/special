#include <cstdlib>
#include <limits>
#include <stdexcept>
#include <utility>
#include <boost/math/special_functions.hpp>
#include "special.hpp"

namespace {
    template<class F>
    auto boost_call(F f) {
        try {
            return f();
        } catch(boost::math::rounding_error&) {
            throw std::domain_error("FIXME: boost::math::rounding_error");
        } catch(boost::math::evaluation_error&) {
            throw std::domain_error("FIXME: boost::math::evaluation_error");
        }
    }
} // unnamed namespace

namespace std {
    double assoc_laguerre(unsigned n, unsigned m, double x) {
        if (std::isnan(x)) return x;
        return boost_call([=]{ return boost::math::laguerre(n, m, x); });
    }

    double assoc_legendre(unsigned l, unsigned m, double x) {
        if (std::isnan(x)) return x;
        return boost_call([=]{ return boost::math::legendre_p(l, m, x); });
    }

    double beta(double x, double y) {
        return boost_call([=]{ return boost::math::beta(x, y); });
    }
    float betaf(float x, float y) {
        return boost_call([=]{ return boost::math::beta(x, y); });
    }

    double legendre(unsigned l, double x) {
        if (std::isnan(x)) return x;
        return boost_call([=]{ return boost::math::legendre_p(l, x); });
    }

    double hypot(double dx, double dy, double dz) _NOEXCEPT
    {
        auto check = [](double value) {
            switch (std::fpclassify(value))
            {
            case FP_NAN:
            case FP_INFINITE:
                std::abort();
                // "fallthrough" but not really
            default:
                return std::abs(value);
            }
        };

        dx = check(dx);

        dy = check(dy);
        if (dx < dy) std::swap(dx, dy);

        dz = check(dz);
        if (dx < dz) std::swap(dx, dz);

        if (dx == 0.0) return 0.0;

        auto frac = [](double const numerator, double const denominator) {
            double const result = numerator / denominator;
            return (result * result);
        };
        return (dx * std::sqrt(frac(dy, dx) + frac(dz, dx) + 1));
    }
} // namespace std
