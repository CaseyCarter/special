#define BOOST_TEST_MODULE beta
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

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
        auto const eps = boost::math::tools::epsilon<T>();

        {
#include "math/test/beta_small_data.ipp"
            for(auto const& datum : beta_small_data) {
                BOOST_CHECK_CLOSE_FRACTION(f(datum[0], datum[1]), datum[2], 3 * eps);
            }
        }

        {
#include "math/test/beta_med_data.ipp"
            for(auto const& datum : beta_med_data) {
                BOOST_CHECK_CLOSE_FRACTION(f(datum[0], datum[1]), datum[2], 98 * eps);
            }
        }

        {
#include "math/test/beta_exp_data.ipp"
            for(auto const& datum : beta_exp_data) {
                BOOST_CHECK_CLOSE_FRACTION(f(datum[0], datum[1]), datum[2], 12 * eps);
            }
        }
#undef SC_
    }
} // unnamed namespace

BOOST_AUTO_TEST_CASE(beta)  { test_beta<double>(std::beta); }
BOOST_AUTO_TEST_CASE(betaf) { test_beta<float>(std::betaf); }
BOOST_AUTO_TEST_CASE(betal) { test_beta<long double>(std::betal); }

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
