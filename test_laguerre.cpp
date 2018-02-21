#define BOOST_TEST_MODULE laguerre
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
    void test_laguerre(F f) {
#define SC_(X) static_cast<T>(BOOST_JOIN(X, L))
        auto const eps = boost::math::tools::epsilon<T>();

        {
#include "math/test/laguerre2.ipp"
            for(auto const& datum : laguerre2) {
                BOOST_CHECK_CLOSE_FRACTION(f(std::lround(datum[0]), datum[1]), datum[2], 4000 * eps);
            }
        }

#if 0
        {
#include "math/test/laguerre3.ipp"
            for(auto const& datum : beta_med_data) {
                BOOST_CHECK_CLOSE_FRACTION(f(datum[0], datum[1]), datum[2], 98 * eps);
            }
        }
#endif

#undef SC_
    }
} // unnamed namespace

BOOST_AUTO_TEST_CASE(laguerre)  { test_laguerre<double>(std::laguerre); }
BOOST_AUTO_TEST_CASE(laguerref) { test_laguerre<float>(std::laguerref); }
BOOST_AUTO_TEST_CASE(laguerrel) { test_laguerre<long double>(std::laguerrel); }

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
