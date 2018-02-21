#define BOOST_TEST_MODULE beta

#include <cmath>
#include <limits>
#include <boost/array.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/test/included/unit_test.hpp>
#include "special.hpp"

namespace {
    template<class T>
    struct table_type { using type = T; };

    template<class T, class F>
    void test_beta(F f) {
#define SC_(X) static_cast<T>(BOOST_JOIN(X, L))
        auto const eps = boost::math::tools::epsilon<double>();

        {
#include "math/test/beta_small_data.ipp"
            for(auto const& datum : beta_small_data) {
                BOOST_CHECK_CLOSE(f(datum[0], datum[1]), datum[2], 3 * eps * 100);
            }
        }

        {
#include "math/test/beta_med_data.ipp"
            for(auto const& datum : beta_med_data) {
                BOOST_CHECK_CLOSE(f(datum[0], datum[1]), datum[2], 98 * eps * 100);
            }
        }

        {
#include "math/test/beta_exp_data.ipp"
            for(auto const& datum : beta_exp_data) {
                BOOST_CHECK_CLOSE(f(datum[0], datum[1]), datum[2], 12 * eps * 100);
            }
        }
#undef SC_
    }
} // unnamed namespace

BOOST_AUTO_TEST_CASE(beta) {
    test_beta<double>([](auto... args){ return std::beta(args...); });
}

BOOST_AUTO_TEST_CASE(betaf) {
    test_beta<float>([](auto... args){ return std::betaf(args...); });
}

BOOST_AUTO_TEST_CASE(betal) {
    test_beta<long double>([](auto... args){ return std::betal(args...); });
}
