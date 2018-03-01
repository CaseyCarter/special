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
    template<class T>
    void single_check(T const x, T const y, T const z, T const result) {
        if (std::isnan(result)) {
            BOOST_CHECK(std::isnan(std::hypot(x, y, z)));
        } else {
            BOOST_CHECK_EQUAL(std::hypot(x, y, z), result);
        }
        BOOST_CHECK(verify_not_domain_error());
    }

    template<class T>
    void permute(T const x, T const y, T const z, T const result) {
        single_check(x, y, z, result);
        single_check(x, z, y, result);
        single_check(y, x, z, result);
        single_check(y, z, x, result);
        single_check(z, x, y, result);
        single_check(z, y, x, result);
    };

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_hypot, T, fptypes) {
        single_check(T{0}, T{0}, T{0}, T{0});
        permute(T{1}, T{0}, T{0}, T{1});
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_hypot_boundaries, T, fptypes) {
        errno = 0;
         // C11 F.10.4.3: "hypot(+/-inf, y) returns +inf even if y is NaN"
        permute(+inf<T>,    T{0}, T{1}, inf<T>);
        permute(-inf<T>,    T{0}, T{1}, inf<T>);
        permute(+inf<T>, qNaN<T>, T{1}, inf<T>);
        permute(-inf<T>, qNaN<T>, T{1}, inf<T>);

        // NaN with no infinity produces NaN
        permute(qNaN<T>, T{0}, T{0}, qNaN<T>);
    }
} // namespace hypot_

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
