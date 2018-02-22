#define BOOST_TEST_MODULE special
#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API

#include <cerrno>
#include <cmath>
#include <limits>
#include <utility>
#include <boost/array.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/precision.hpp>
#include <boost/test/included/unit_test.hpp>
#include "special.hpp"

#undef small // Thanks for this, Windows SDK.
#define SC_(X) static_cast<T>(BOOST_JOIN(X, L))

template<class T>
constexpr auto eps = boost::math::tools::epsilon<T>();

template<class T>
constexpr auto qNaN = std::numeric_limits<T>::quiet_NaN();

template<class T>
struct table_type { using type = T; };

template<class>
constexpr bool always_false = false;

using fptypes = boost::mpl::list<float, double, long double>;

bool verify_domain_error() {
    return std::exchange(errno, 0) == EDOM;
}

bool verify_not_domain_error() {
    return std::exchange(errno, 0) == 0;
}

namespace assoc_laguerre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>, "BOOM");
    };
    template<>
    constexpr auto test_fn<float> = std::assoc_laguerref;
    template<>
    constexpr auto test_fn<double> = std::assoc_laguerre;
    template<>
    constexpr auto test_fn<long double> = std::assoc_laguerrel;

    template<class T>
    constexpr auto control_fn = [](unsigned n, unsigned m, T x) {
        return boost::math::laguerre(n, m, x);
    };

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_laguerre, T, fptypes) {
#include "math/test/laguerre3.ipp"

        for(auto const& datum : laguerre3) {
            unsigned const n = std::lround(datum[0]);
            unsigned const m = std::lround(datum[1]);
            auto const actual = test_fn<T>(n, m, datum[2]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(n, m, datum[2]));
            if (actual != datum[3]) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[3], 440 * eps<T>);
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_laguerre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, 5, static_cast<T>(0.5L)), static_cast<T>(88.31510416666666666666666666666666666667L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 0, static_cast<T>(2.5L)), static_cast<T>(-0.8802526766660982969576719576719576719577L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 1, static_cast<T>(4.5L)), static_cast<T>(1.564311458042689732142857142857142857143L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 6, static_cast<T>(8.5L)), static_cast<T>(20.51596541066649098875661375661375661376L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 12, static_cast<T>(12.5L)), static_cast<T>(-199.5560968456234671241181657848324514991L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, 40, static_cast<T>(12.5L)), static_cast<T>(-4.996769495006119488583146995907246595400e16L), tolerance);
    }
} // namespace assoc_laguerre

namespace assoc_legendre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>, "BOOM");
    };
    template<>
    constexpr auto test_fn<float> = std::assoc_legendref;
    template<>
    constexpr auto test_fn<double> = std::assoc_legendre;
    template<>
    constexpr auto test_fn<long double> = std::assoc_legendrel;

    template<class T>
    constexpr auto control_fn = [](unsigned l, unsigned m, T x) {
        return boost::math::legendre_p(l, m, x);
    };

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_legendre, T, fptypes) {
#include "math/test/assoc_legendre_p.ipp"

        for(auto const& datum : assoc_legendre_p) {
            auto const l = static_cast<unsigned const>(datum[0]);
            unsigned const m = std::lround(datum[1]);
            auto const actual = test_fn<T>(l, m, datum[2]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(l, m, datum[2]));
            if (actual != datum[3]) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[3], 150 * eps<T>); // FIXME: tune
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_legendre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, 2, static_cast<T>(0.5L)), static_cast<T>(4.218750000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, 2, static_cast<T>(0.5L)), static_cast<T>(5.625000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, 5, static_cast<T>(0.5L)), static_cast<T>(-5696.789530152175143607977274672800795328L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, 4, static_cast<T>(0.5L)), static_cast<T>(465.1171875000000000000000000000000000000L), tolerance);
        if(std::numeric_limits<T>::max_exponent > std::numeric_limits<float>::max_exponent) {
            BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, 30, static_cast<T>(0.5L)), static_cast<T>(-7.855722083232252643913331343916012143461e45L), tolerance);
        }
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, 20, static_cast<T>(0.5L)), static_cast<T>(4.966634149702370788037088925152355134665e30L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, 2, static_cast<T>(-0.5L)), static_cast<T>(4.218750000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, 2, static_cast<T>(-0.5L)), static_cast<T>(-5.625000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, 5, static_cast<T>(-0.5L)), static_cast<T>(-5696.789530152175143607977274672800795328L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, 4, static_cast<T>(-0.5L)), static_cast<T>(465.1171875000000000000000000000000000000L), tolerance);
        if(std::numeric_limits<T>::max_exponent > std::numeric_limits<float>::max_exponent) {
            BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, 30, static_cast<T>(-0.5L)), static_cast<T>(-7.855722083232252643913331343916012143461e45L), tolerance);
        }
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, 20, static_cast<T>(-0.5L)), static_cast<T>(-4.966634149702370788037088925152355134665e30L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, -2, static_cast<T>(0.5L)), static_cast<T>(0.01171875000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, -2, static_cast<T>(0.5L)), static_cast<T>(0.04687500000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, -5, static_cast<T>(0.5L)), static_cast<T>(0.00002378609812640364935569308025139290054701L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, -4, static_cast<T>(0.5L)), static_cast<T>(0.0002563476562500000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, -30, static_cast<T>(0.5L)), static_cast<T>(-2.379819988646847616996471299410611801239e-48L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, -20, static_cast<T>(0.5L)), static_cast<T>(4.356454600748202401657099008867502679122e-33L), tolerance);

        BOOST_CHECK_CLOSE_FRACTION(::boost::math::legendre_q(1, static_cast<T>(0.5L)), static_cast<T>(-0.7253469278329725771511886907693685738381L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(::boost::math::legendre_q(4, static_cast<T>(0.5L)), static_cast<T>(0.4401745259867706044988642951843745400835L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(::boost::math::legendre_q(7, static_cast<T>(0.5L)), static_cast<T>(-0.3439152932669753451878700644212067616780L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(::boost::math::legendre_q(40, static_cast<T>(0.5L)), static_cast<T>(0.1493671665503550095010454949479907886011L), tolerance);
    }
} // namespace assoc_legendre

namespace beta {
    template<class T>
    constexpr auto control_fn = [](T x, T y) {
        return boost::math::beta(x, y);
    };

    template<class T>
    constexpr auto test_fn = [](auto x, T) {
        static_assert(always_false<decltype(x)>, "BOOM");
    };
    template<>
    constexpr auto test_fn<float> = std::betaf;
    template<>
    constexpr auto test_fn<double> = std::beta;
    template<>
    constexpr auto test_fn<long double> = std::betal;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_beta, T, fptypes) {
#include "math/test/beta_small_data.ipp"
#include "math/test/beta_med_data.ipp"
#include "math/test/beta_exp_data.ipp"

        for(auto const& datum : beta_small_data) {
            auto const actual = test_fn<T>(datum[0], datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
            BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 3 * eps<T>);
        }

        for(auto const& datum : beta_med_data) {
            auto const actual = test_fn<T>(datum[0], datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
            BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 98 * eps<T>);
        }

        for(auto const& datum : beta_exp_data) {
            auto const actual = test_fn<T>(datum[0], datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
            BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 12 * eps<T>);
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_beta_spots, T, fptypes) {
        auto const tolerance = eps<T> * 20;
        auto const small = eps<T> / 1024;

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1), static_cast<T>(1)), static_cast<T>(1), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1), static_cast<T>(4)), static_cast<T>(0.25), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(4), static_cast<T>(1)), static_cast<T>(0.25), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(small, static_cast<T>(4)), 1/small, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(4), small), 1/small, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(4), static_cast<T>(20)), static_cast<T>(0.00002823263692828910220214568040654997176736L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(0.0125L), static_cast<T>(0.000023L)), static_cast<T>(43558.24045647538375006349016083320744662L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_beta_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(0), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace beta

namespace laguerre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>, "BOOM");
    };
    template<>
    constexpr auto test_fn<float> = std::laguerref;
    template<>
    constexpr auto test_fn<double> = std::laguerre;
    template<>
    constexpr auto test_fn<long double> = std::laguerrel;

    template<class T>
    constexpr auto control_fn = [](unsigned n, T x) {
        return boost::math::laguerre(n, x);
    };

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_laguerre, T, fptypes) {
#include "math/test/laguerre2.ipp"

        for(auto const& datum : laguerre2) {
            unsigned const n = std::lround(datum[0]);
            auto const actual = test_fn<T>(n, datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(n, datum[1]));
            if (actual != datum[2]) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 3200 * eps<T>);
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_laguerre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(0.5L)), static_cast<T>(0.5L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(0.5L)), static_cast<T>(-0.3307291666666666666666666666666666666667L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(0.5L)), static_cast<T>(-0.5183392237103174603174603174603174603175L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(20, static_cast<T>(0.5L)), static_cast<T>(0.3120174870800154148915399248893113634676L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, static_cast<T>(0.5L)), static_cast<T>(-0.3181388060269979064951118308575628226834L), tolerance);

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(-0.5L)), static_cast<T>(1.5L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(-0.5L)), static_cast<T>(3.835937500000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(-0.5L)), static_cast<T>(7.950934709821428571428571428571428571429L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(20, static_cast<T>(-0.5L)), static_cast<T>(76.12915699869631476833699787070874048223L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, static_cast<T>(-0.5L)), static_cast<T>(2307.428631277506570629232863491518399720L), tolerance);

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(4.5L)), static_cast<T>(-3.500000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(4.5L)), static_cast<T>(0.08593750000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(4.5L)), static_cast<T>(-1.036928013392857142857142857142857142857L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(20, static_cast<T>(4.5L)), static_cast<T>(1.437239150257817378525582974722170737587L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, static_cast<T>(4.5L)), static_cast<T>(-0.7795068145562651416494321484050019245248L), tolerance);
    }
} // namespace laguerre

namespace legendre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>, "BOOM");
    };
    template<>
    constexpr auto test_fn<float> = std::legendref;
    template<>
    constexpr auto test_fn<double> = std::legendre;
    template<>
    constexpr auto test_fn<long double> = std::legendrel;

    template<class T>
    constexpr auto control_fn = [](unsigned l, T x) {
        return boost::math::legendre_p(l, x);
    };

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_legendre, T, fptypes) {
#include "math/test/legendre_p.ipp"
#include "math/test/legendre_p_large.ipp"

        auto test_one_value = [](auto const& datum) {
            unsigned const l = std::lround(datum[0]);
            auto const actual = test_fn<T>(l, datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(l, datum[1]));
            if (actual != datum[2]) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 310 * eps<T>);
        };

        std::for_each(std::begin(legendre_p), std::end(legendre_p), test_one_value);
        std::for_each(std::begin(legendre_p_large), std::end(legendre_p_large), test_one_value);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_legendre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(0.5L)), static_cast<T>(0.5L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-1, static_cast<T>(0.5L)), static_cast<T>(1L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(0.5L)), static_cast<T>(-0.2890625000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, static_cast<T>(0.5L)), static_cast<T>(-0.4375000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(0.5L)), static_cast<T>(0.2231445312500000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, static_cast<T>(0.5L)), static_cast<T>(0.3232421875000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, static_cast<T>(0.5L)), static_cast<T>(-0.09542943523261546936538467572384923220258L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, static_cast<T>(0.5L)), static_cast<T>(-0.1316993126940266257030910566308990611306L), tolerance);
    }
} // namespace legendre

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
