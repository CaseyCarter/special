#define BOOST_TEST_MODULE legendre
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <cmath>
#include <algorithm>
#include <limits>
#include <boost/array.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/test/included/unit_test.hpp>
#include "special.hpp"

BOOST_AUTO_TEST_CASE(nothing) {}

template<class T>
constexpr auto eps = boost::math::tools::epsilon<T>();

template<class T>
struct table_type { using type = T; };

using fptypes = boost::mpl::list<float, double, long double>;

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
