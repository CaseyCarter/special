#define BOOST_TEST_MODULE beta
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
    constexpr bool always_false = false;

    template<class T>
    constexpr auto control_fn = [](T x, T y) { return boost::math::beta(x, y); };

    template<class T>
    constexpr auto test_fn = [](auto x, T) { static_assert(always_false<decltype(x)>, "BOOM"); };
    template<>
    constexpr auto test_fn<float> = std::betaf;
    template<>
    constexpr auto test_fn<double> = std::beta;
    template<>
    constexpr auto test_fn<long double> = std::betal;

    using types = boost::mpl::list<float, double, long double>;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_beta, T, types) {
#define SC_(X) static_cast<T>(BOOST_JOIN(X, L))
#include "math/test/beta_small_data.ipp"
#include "math/test/beta_med_data.ipp"
#include "math/test/beta_exp_data.ipp"
#undef SC_

        auto const eps = boost::math::tools::epsilon<T>();

        for(auto const& datum : beta_small_data) {
            auto const actual = test_fn<T>(datum[0], datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
            BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 3 * eps);
        }

        for(auto const& datum : beta_med_data) {
            auto const actual = test_fn<T>(datum[0], datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
            BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 98 * eps);
        }

        for(auto const& datum : beta_exp_data) {
            auto const actual = test_fn<T>(datum[0], datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
            BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 12 * eps);
        }
    }
} // unnamed namespace

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
