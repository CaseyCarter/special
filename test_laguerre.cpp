#define BOOST_TEST_MODULE laguerre
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <cmath>
#include <limits>
#include <boost/array.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/test/included/unit_test.hpp>
#include "special.hpp"

namespace {
    template<class T>
    struct table_type { using type = T; };

    template<class>
    constexpr auto test_fn = [](unsigned int, auto x) { static_assert(always_false<decltype(x)>, "BOOM"); };
    template<>
    constexpr auto test_fn<float> = std::laguerref;
    template<>
    constexpr auto test_fn<double> = std::laguerre;
    template<>
    constexpr auto test_fn<long double> = std::laguerrel;

    template<class T>
    constexpr auto control_fn = [](unsigned int n, T x) { return boost::math::laguerre(n, x); };

    using types = boost::mpl::list<float, double, long double>;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_laguerre, T, types) {
#define SC_(X) static_cast<T>(BOOST_JOIN(X, L))
#include "math/test/laguerre2.ipp"
#include "math/test/laguerre3.ipp"
#undef SC_

        auto const eps = boost::math::tools::epsilon<T>();

        for(auto const& datum : laguerre2) {
            unsigned int const n = std::lround(datum[0]);
            auto const actual = test_fn<T>(n, datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(n, datum[1]));
            BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 4000 * eps);
        }

#if 0
        for(auto const& datum : beta_med_data) {
            BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(datum[0], datum[1]), datum[2], 98 * eps);
        }
#endif
    }
} // unnamed namespace

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
