#define BOOST_TEST_MODULE WIP
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <cmath>
#include <algorithm>
#include <array>
#include <limits>
#include <boost/array.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/test/included/unit_test.hpp>
#include "special.hpp"

BOOST_AUTO_TEST_CASE(nothing) {}

#undef small // Thanks for this, Windows SDK.
#define SC_(X) static_cast<T>(BOOST_JOIN(X, L))

template<class T>
constexpr auto eps = boost::math::tools::epsilon<T>();

template<class T>
constexpr auto qNaN = std::numeric_limits<T>::quiet_NaN();

template<class T>
constexpr auto inf = std::numeric_limits<T>::infinity();

template<class T>
struct table_type { using type = T; };

template<class>
constexpr bool always_false = false;

using fptypes = boost::mpl::list<float, double, long double>;

inline bool verify_domain_error() {
    return std::exchange(errno, 0) == EDOM;
}

inline bool verify_not_domain_error() {
    return std::exchange(errno, 0) == 0;
}

template<class Range, class Fn>
inline void for_each(Range&& rng, Fn f) {
    using std::begin;
    using std::end;
    std::for_each(begin(rng), end(rng), f);
}

template<class Range1, class Range2>
inline bool equal(Range1&& r1, Range2&& r2) {
    using std::begin;
    using std::end;
    return std::equal(begin(r1), end(r1), begin(r2), end(r2));
}

namespace hypot_ {
    BOOST_AUTO_TEST_CASE_TEMPLATE(test_hypot, T, fptypes) {
        auto test_hypot = [](auto x, auto y, auto z, auto result) {
            BOOST_CHECK_EQUAL(std::hypot(x, y, z), result);
            BOOST_CHECK_EQUAL(std::hypot(x, z, y), result);
            BOOST_CHECK_EQUAL(std::hypot(y, x, z), result);
            BOOST_CHECK_EQUAL(std::hypot(y, z, x), result);
            BOOST_CHECK_EQUAL(std::hypot(z, x, y), result);
            BOOST_CHECK_EQUAL(std::hypot(z, y, x), result);
        };

        test_hypot(0.0, 0.0, 0.0, 0.0);
        test_hypot(1.0, 0.0, 0.0, 1.0);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_hypot_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(std::hypot(qNaN<T>, 0.0, 0.0)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(std::hypot(0.0, qNaN<T>, 0.0)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(std::hypot(0.0, 0.0, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        auto test_hypot = [](auto x, auto y, auto z, auto result) {
            BOOST_CHECK_EQUAL(std::hypot(x, y, z), result);
            BOOST_CHECK(verify_not_domain_error());
            BOOST_CHECK_EQUAL(std::hypot(x, z, y), result);
            BOOST_CHECK(verify_not_domain_error());
            BOOST_CHECK_EQUAL(std::hypot(y, x, z), result);
            BOOST_CHECK(verify_not_domain_error());
            BOOST_CHECK_EQUAL(std::hypot(y, z, x), result);
            BOOST_CHECK(verify_not_domain_error());
            BOOST_CHECK_EQUAL(std::hypot(z, x, y), result);
            BOOST_CHECK(verify_not_domain_error());
            BOOST_CHECK_EQUAL(std::hypot(z, y, x), result);
            BOOST_CHECK(verify_not_domain_error());
        };

        test_hypot(+inf<T>,     0.0, 1.0, inf<T>);
        test_hypot(-inf<T>,     0.0, 1.0, inf<T>);
        test_hypot(+inf<T>, qNaN<T>, 1.0, inf<T>); // C11 F.10.4.3: "hypot(+/-inf, y)
        test_hypot(-inf<T>, qNaN<T>, 1.0, inf<T>); // returns +inf even if y is NaN"
    }
} // namespace hypot_

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
