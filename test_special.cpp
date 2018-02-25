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

namespace assoc_laguerre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>);
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
            if (!(actual == datum[3])) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[3], 440 * eps<T>);
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_laguerre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, 5, static_cast<T>(0.5L)),
            static_cast<T>(88.31510416666666666666666666666666666667L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 0, static_cast<T>(2.5L)),
            static_cast<T>(-0.8802526766660982969576719576719576719577L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 1, static_cast<T>(4.5L)),
            static_cast<T>(1.564311458042689732142857142857142857143L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 6, static_cast<T>(8.5L)),
            static_cast<T>(20.51596541066649098875661375661375661376L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, 12, static_cast<T>(12.5L)),
            static_cast<T>(-199.5560968456234671241181657848324514991L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, 40, static_cast<T>(12.5L)),
            static_cast<T>(-4.996769495006119488583146995907246595400e16L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_laguerre_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(0u, 0u, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace assoc_laguerre

namespace assoc_legendre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>);
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
            if (!(actual == datum[3])) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[3], 150 * eps<T>); // FIXME: tune
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_legendre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, 2, static_cast<T>(0.5L)),
            static_cast<T>(4.218750000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, 2, static_cast<T>(0.5L)),
            static_cast<T>(5.625000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, 5, static_cast<T>(0.5L)),
            static_cast<T>(-5696.789530152175143607977274672800795328L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, 4, static_cast<T>(0.5L)),
            static_cast<T>(465.1171875000000000000000000000000000000L), tolerance);
        if(std::numeric_limits<T>::max_exponent > std::numeric_limits<float>::max_exponent) {
            BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, 30, static_cast<T>(0.5L)),
                static_cast<T>(-7.855722083232252643913331343916012143461e45L), tolerance);
        }
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, 20, static_cast<T>(0.5L)),
            static_cast<T>(4.966634149702370788037088925152355134665e30L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, 2, static_cast<T>(-0.5L)),
            static_cast<T>(4.218750000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, 2, static_cast<T>(-0.5L)),
            static_cast<T>(-5.625000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, 5, static_cast<T>(-0.5L)),
            static_cast<T>(-5696.789530152175143607977274672800795328L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, 4, static_cast<T>(-0.5L)),
            static_cast<T>(465.1171875000000000000000000000000000000L), tolerance);
        if(std::numeric_limits<T>::max_exponent > std::numeric_limits<float>::max_exponent) {
            BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, 30, static_cast<T>(-0.5L)),
                static_cast<T>(-7.855722083232252643913331343916012143461e45L), tolerance);
        }
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, 20, static_cast<T>(-0.5L)),
            static_cast<T>(-4.966634149702370788037088925152355134665e30L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, -2, static_cast<T>(0.5L)),
            static_cast<T>(0.01171875000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, -2, static_cast<T>(0.5L)),
            static_cast<T>(0.04687500000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, -5, static_cast<T>(0.5L)),
            static_cast<T>(0.00002378609812640364935569308025139290054701L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, -4, static_cast<T>(0.5L)),
            static_cast<T>(0.0002563476562500000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, -30, static_cast<T>(0.5L)),
            static_cast<T>(-2.379819988646847616996471299410611801239e-48L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, -20, static_cast<T>(0.5L)),
            static_cast<T>(4.356454600748202401657099008867502679122e-33L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_assoc_legendre_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(4, 2, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |x| <= 1
        BOOST_CHECK_EQUAL(test_fn<T>(4, 2, static_cast<T>(1)), static_cast<T>(0));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_EQUAL(test_fn<T>(4, 2, static_cast<T>(-1)), static_cast<T>(0));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(4, 2, static_cast<T>(32))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(4, 2, static_cast<T>(-32))));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace assoc_legendre

namespace beta {
    template<class T>
    constexpr auto control_fn = [](T x, T y) {
        return boost::math::beta(x, y);
    };

    template<class T>
    constexpr auto test_fn = [](auto x, T) {
        static_assert(always_false<decltype(x)>);
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

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0], datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
            };
        };

        ::for_each(beta_small_data, tester(3 * eps<T>));
        ::for_each(beta_med_data, tester(98 * eps<T>));
        ::for_each(beta_exp_data, tester(12 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_beta_spots, T, fptypes) {
        auto const tolerance = eps<T> * 20;
        auto const small = eps<T> / 1024;

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1), static_cast<T>(1)),
            static_cast<T>(1), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1), static_cast<T>(4)),
            static_cast<T>(0.25), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(4), static_cast<T>(1)),
            static_cast<T>(0.25), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(small, static_cast<T>(4)), 1/small, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(4), small), 1/small, tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(4), static_cast<T>(20)),
            static_cast<T>(0.00002823263692828910220214568040654997176736L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(0.0125L), static_cast<T>(0.000023L)),
            static_cast<T>(43558.24045647538375006349016083320744662L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_beta_boundaries, T, fptypes) {
        auto const tolerance = eps<T> * 20;

        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is x > 0 && y > 0
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(0), static_cast<T>(1))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(0), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1e-20L), static_cast<T>(1e-20L)),
            static_cast<T>(1.99999999999999999999999999999999999999967101e20L), tolerance);
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace beta

namespace comp_ellint_1 {
    template<class T>
    constexpr auto control_fn = [](T k) {
        return boost::math::ellint_1(k);
    };

    template<class T>
    constexpr auto test_fn = [](auto k) {
        static_assert(always_false<decltype(k)>);
    };
    template<>
    constexpr auto test_fn<float> = std::comp_ellint_1f;
    template<>
    constexpr auto test_fn<double> = std::comp_ellint_1;
    template<>
    constexpr auto test_fn<long double> = std::comp_ellint_1l;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_comp_ellint_1, T, fptypes) {
        // From function test_spots in test_ellint_1.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 2>, 9> data2 = {{
            {{ SC_(0.0), SC_(1.5707963267948966192313216916397514420985846996876) }},
            {{ SC_(0.125), SC_(1.5769867712158131421244030532288080803822271060839) }},
            {{ SC_(0.25), SC_(1.5962422221317835101489690714979498795055744578951) }},
            {{ SC_(0.29296875) /*T(300)/1024*/, SC_(1.6062331054696636704261124078746600894998873503208) }},
            {{ SC_(0.390625) /*T(400)/1024*/, SC_(1.6364782007562008756208066125715722889067992997614) }},
            {{ SC_(-0.5), SC_(1.6857503548125960428712036577990769895008008941411) }},
            {{ SC_(-0.75), SC_(1.9109897807518291965531482187613425592531451316788) }},
            {{ SC_(0.875) /*1-T(1)/8*/, SC_(2.185488469278223686913080323730158689730428415766) }},
            {{ SC_(0.9990234375) /*1-T(1)/1024*/, SC_(4.5074135978990422666372495313621124487894807327687) }},
        }};

#include "math/test/ellint_k_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0]));
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[1], tolerance);
            };
        };

        ::for_each(data2, tester(eps<T>));
        ::for_each(ellint_k_data, tester(eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_comp_ellint_1_boundaries, T, fptypes) {
        auto const tolerance = eps<T>;

        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |k| <= 1
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(-2))));
        BOOST_CHECK(verify_domain_error());

        BOOST_CHECK_EQUAL(test_fn<T>(static_cast<T>(1)), inf<T>);
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_EQUAL(test_fn<T>(static_cast<T>(-1)), inf<T>);
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace comp_ellint_1

namespace comp_ellint_2 {
    template<class T>
    constexpr auto control_fn = [](T k) {
        return boost::math::ellint_2(k);
    };

    template<class T>
    constexpr auto test_fn = [](auto k) {
        static_assert(always_false<decltype(k)>);
    };
    template<>
    constexpr auto test_fn<float> = std::comp_ellint_2f;
    template<>
    constexpr auto test_fn<double> = std::comp_ellint_2;
    template<>
    constexpr auto test_fn<long double> = std::comp_ellint_2l;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_comp_ellint_2, T, fptypes) {
        // From function test_spots in test_ellint_2.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 2>, 10> data2 = {{
            {{ SC_(-1.0), SC_(1.0) }},
            {{ SC_(0.0), SC_(1.5707963267948966192313216916397514420985846996876) }},
            {{ SC_(0.09765625) /*T(100) / 1024*/, SC_(1.5670445330545086723323795143598956428788609133377) }},
            {{ SC_(0.1953125) /*T(200) / 1024*/, SC_(1.5557071588766556854463404816624361127847775545087) }},
            {{ SC_(0.29296875) /*T(300) / 1024*/, SC_(1.5365278991162754883035625322482669608948678755743) }},
            {{ SC_(0.390625) /*T(400) / 1024*/, SC_(1.5090417763083482272165682786143770446401437564021) }},
            {{ SC_(-0.5), SC_(1.4674622093394271554597952669909161360253617523272) }},
            {{ SC_(-0.5859375) /*T(-600) / 1024*/, SC_(1.4257538571071297192428217218834579920545946473778) }},
            {{ SC_(-0.78125) /*T(-800) / 1024*/, SC_(1.2927868476159125056958680222998765985004489572909) }},
            {{ SC_(-0.87890625) /*T(-900) / 1024*/, SC_(1.1966864890248739524112920627353824133420353430982) }},
        }};

#include "math/test/ellint_e_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0]));
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[1], tolerance);
            };
        };

        ::for_each(data2, tester(1.5 * eps<T>));
        ::for_each(ellint_e_data, tester(2 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_comp_ellint_2_boundaries, T, fptypes) {
        auto const tolerance = eps<T>;

        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |k| <= 1
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(-2))));
        BOOST_CHECK(verify_domain_error());

        BOOST_CHECK_EQUAL(test_fn<T>(static_cast<T>(1)), static_cast<T>(1));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_EQUAL(test_fn<T>(static_cast<T>(-1)), static_cast<T>(1));
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace comp_ellint_2

namespace comp_ellint_3 {
    template<class T>
    constexpr auto control_fn = [](T k, T nu) {
        return boost::math::ellint_3(k, nu);
    };

    template<class T>
    constexpr auto test_fn = [](auto k, T) {
        static_assert(always_false<decltype(k)>);
    };
    template<>
    constexpr auto test_fn<float> = std::comp_ellint_3f;
    template<>
    constexpr auto test_fn<double> = std::comp_ellint_3;
    template<>
    constexpr auto test_fn<long double> = std::comp_ellint_3l;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_comp_ellint_3, T, fptypes) {
        // From function test_spots in test_ellint_3.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 17> data2 = {{
            {{ SC_(0.0), SC_(0.2), SC_(1.586867847454166237308008033828114192951) }},
            {{ SC_(0.0), SC_(0.4), SC_(1.639999865864511206865258329748601457626) }},
            {{ SC_(0.0), SC_(0.0), SC_(1.57079632679489661923132169163975144209858469968755291048747) }},
            {{ SC_(0.5), SC_(0.0), SC_(2.221441469079183123507940495030346849307) }},
            {{ SC_(-4.0), SC_(0.3), SC_(0.712708870925620061597924858162260293305195624270730660081949) }},
            {{ SC_(-1e+05), SC_(-0.5), SC_(0.00496944596485066055800109163256108604615568144080386919012831) }},
            {{ SC_(-1e+10), SC_(-0.75), SC_(0.0000157080225184890546939710019277357161497407143903832703317801) }},
            {{ SC_(0.0009765625) /*T(1) / 1024*/, SC_(-0.875), SC_(2.18674503176462374414944618968850352696579451638002110619287) }},
            {{ SC_(0.9990234375) /*T(1023)/1024*/, SC_(-0.875), SC_(101.045289804941384100960063898569538919135722087486350366997) }},
            // Bug cases from Rocco Romeo:
            { { SC_(1e-175), SC_(0.0), SC_(1.57079632679489661923132169163975144209858469968755291048747) } },
            { { SC_(1e-170), SC_(1E-164), SC_(1.57079632679489661923132169163975144209858469968755291048747) } },
            { { SC_(1e-170), SC_(-1E-164), SC_(1.57079632679489661923132169163975144209858469968755291048747) } },
            { { SC_(-3.3306690738754696212708950042724609375e-16) /*-1.5f * ldexp(T(1), -52)*/, SC_(-0.9375), SC_(2.48840049140103464299631535211815755485846563527849342319632) } },
            { { SC_(-3.3306690738754696212708950042724609375e-16) /*-1.5f * ldexp(T(1), -52)*/, SC_(0.9375), SC_(2.48840049140103464299631535211815755485846563527849342319632) } },
            { { SC_(2.6497349136889904540094297102338792048842007654486987513108141219242520439074860067791572114971742036578602e-169) /*ldexp(T(1), -560)*/, SC_(2.1382117680737565169124291737211855030521575040840389583695499937283189127897042869363986028474755585193634e-50) /*ldexp(T(1), -165)*/, SC_(1.57079632679489661923132169163975144209858469968756130722545) } },
            { { SC_(2.6497349136889904540094297102338792048842007654486987513108141219242520439074860067791572114971742036578602e-169) /*ldexp(T(1), -560)*/, SC_(-2.1382117680737565169124291737211855030521575040840389583695499937283189127897042869363986028474755585193634e-50) /*-ldexp(T(1), -165)*/, SC_(1.57079632679489661923132169163975144209858469968754451374949) } },
            { { std::numeric_limits<T>::max_exponent > 600 ? SC_(-4.14951556888099295851240786369116115101244623224243689999565732969065281141290814639970704894710379428819788e180) /*T(-ldexp(T(1), 600))*/ : SC_(0.0), SC_(0.5), std::numeric_limits<T>::max_exponent > 600 ? SC_(7.71118598318249916481121898327895181916104121635240801895419e-91) : SC_(1.68575035481259604287120365779907698950080089414108904411995) } },
        } };

#include "math/test/ellint_pi2_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[1], datum[0]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[1], datum[0]));
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
            };
        };

        ::for_each(data2, tester(eps<T>));
        ::for_each(ellint_pi2_data, tester(3 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_comp_ellint_3_boundaries, T, fptypes) {
        auto const tolerance = eps<T>;

        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(0), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(0.75))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |k| < 1
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), static_cast<T>(0.75))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(-1), static_cast<T>(0.75))));
        BOOST_CHECK(verify_domain_error());

        BOOST_CHECK(!std::isnan(test_fn<T>(static_cast<T>(1) - eps<T>, static_cast<T>(0.75))));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(!std::isnan(test_fn<T>(static_cast<T>(-1) + eps<T>, static_cast<T>(0.75))));
        BOOST_CHECK(verify_not_domain_error());

        BOOST_CHECK(std::isnan(test_fn<T>(T(1.0001), T(-1))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(0.5), T(1))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(0.5), T(2))));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace comp_ellint_3

namespace cyl_bessel_i {
    template<class T>
    constexpr auto control_fn = [](T nu, T x) {
        return boost::math::cyl_bessel_i(nu, x);
    };

    template<class T>
    constexpr auto test_fn = [](auto nu, T) {
        static_assert(always_false<decltype(nu)>);
    };
    template<>
    constexpr auto test_fn<float> = std::cyl_bessel_if;
    template<>
    constexpr auto test_fn<double> = std::cyl_bessel_i;
    template<>
    constexpr auto test_fn<long double> = std::cyl_bessel_il;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_i, T, fptypes) {
        // Data taken from test_bessel_i.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 10> i0_data = {{
            {{ SC_(0.0), SC_(0.0), SC_(1.0) }},
            {{ SC_(0.0), SC_(1.0), SC_(1.26606587775200833559824462521471753760767031135496220680814) }},
            {{ SC_(0.0), SC_(-2.0), SC_(2.27958530233606726743720444081153335328584110278545905407084) }},
            {{ SC_(0.0), SC_(4.0), SC_(11.3019219521363304963562701832171024974126165944353377060065) }},
            {{ SC_(0.0), SC_(-7.0), SC_(168.593908510289698857326627187500840376522679234531714193194) }},
            {{ SC_(0.0), SC_(0.0009765625), SC_(1.00000023841859331241759166109699567801556273303717896447683) }},
            {{ SC_(0.0), SC_(9.5367431640625e-7), SC_(1.00000000000022737367544324498417583090700894607432256476338) }},
            {{ SC_(0.0), SC_(-1.0), SC_(1.26606587775200833559824462521471753760767031135496220680814) }},
            {{ SC_(0.0), SC_(100.0), SC_(1.07375170713107382351972085760349466128840319332527279540154e42) }},
            {{ SC_(0.0), SC_(200.0), SC_(2.03968717340972461954167312677945962233267573614834337894328e85) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 10> i1_data = {{
            {{ SC_(1.0), SC_(0.0), SC_(0.0) }},
            {{ SC_(1.0), SC_(1.0), SC_(0.565159103992485027207696027609863307328899621621092009480294) }},
            {{ SC_(1.0), SC_(-2.0), SC_(-1.59063685463732906338225442499966624795447815949553664713229) }},
            {{ SC_(1.0), SC_(4.0), SC_(9.75946515370444990947519256731268090005597033325296730692753) }},
            {{ SC_(1.0), SC_(-8.0), SC_(-399.873136782560098219083086145822754889628443904067647306574) }},
            {{ SC_(1.0), SC_(0.0009765625), SC_(0.000488281308207663226432087816784315537514225208473395063575150) }},
            {{ SC_(1.0), SC_(9.5367431640625e-7), SC_(4.76837158203179210108624277276025646653133998635956784292029E-7) }},
            {{ SC_(1.0), SC_(-10.0), SC_(-2670.98830370125465434103196677215254914574515378753771310849) }},
            {{ SC_(1.0), SC_(100.0), SC_(1.06836939033816248120614576322429526544612284405623226965918e42) }},
            {{ SC_(1.0), SC_(200.0), SC_(2.03458154933206270342742797713906950389661161681122964159220e85) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 11> in_data = {{
            {{ SC_(-2.0), SC_(0.0), SC_(0.0) }},
            {{ SC_(2.0), SC_(9.5367431640625e-7), SC_(1.13686837721624646204093977095674566928522671779753217215467e-13) }},
            {{ SC_(5.0), SC_(10.0), SC_(777.188286403259959907293484802339632852674154572666041953297) }},
            {{ SC_(-5.0), SC_(100.0), SC_(9.47009387303558124618275555002161742321578485033007130107740e41) }},
            {{ SC_(-5.0), SC_(-1.0), SC_(-0.000271463155956971875181073905153777342383564426758143634974124) }},
            {{ SC_(10.0), SC_(20.0), SC_(3.54020020901952109905289138244985607057267103782948493874391e6) }},
            {{ SC_(10.0), SC_(-5.0), SC_(0.00458004441917605126118647027872016953192323139337073320016447) }},
            {{ SC_(1e+02), SC_(9.0), SC_(2.74306601746058997093587654668959071522869282506446891736820e-93) }},
            {{ SC_(1e+02), SC_(80.0), SC_(4.65194832850610205318128191404145885093970505338730540776711e8) }},
            {{ SC_(-100.0), SC_(-200.0), SC_(4.35275044972702191438729017441198257508190719030765213981307e74) }},
            {{ SC_(10.0), SC_(1e-100), SC_(2.69114445546737213403880070546737213403880070546737213403880e-1010) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 10> iv_data = {{
            {{ SC_(2.25), SC_(9.5367431640625e-7), SC_(2.34379212133481347189068464680335815256364262507955635911656e-15) }},
            {{ SC_(5.5), SC_(3.125), SC_(0.0583514045989371500460946536220735787163510569634133670181210) }},
            {{ SC_(-4.9990234375), SC_(2.125), SC_(0.0267920938009571023702933210070984416052633027166975342895062) }},
            {{ SC_(-5.5), SC_(10.0), SC_(597.577606961369169607937419869926705730305175364662688426534) }},
            {{ SC_(-5.5), SC_(100.0), SC_(9.22362906144706871737354069133813819358704200689067071415379e41) }},
            {{ SC_(-10.0002994537353515625), SC_(0.0009765625), SC_(1.41474005665181350367684623930576333542989766867888186478185e35) }},
            {{ SC_(-10.0002994537353515625), SC_(50.0), SC_(1.07153277202900671531087024688681954238311679648319534644743e20) }},
            {{ SC_(141.400390625), SC_(100.0), SC_(2066.27694757392660413922181531984160871678224178890247540320) }},
            {{ SC_(141.400390625), SC_(200.0), SC_(2.23699739472246928794922868978337381373643889659337595319774e64) }},
            {{ SC_(-141.400390625), SC_(100.0), SC_(2066.27694672763190927440969155740243346136463461655104698748) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 5> iv_large_data = {{
            // Bug report https://svn.boost.org/trac/boost/ticket/5560:
            {{ SC_(-1.0), SC_(3.7291703656001033716454826577314669186882357673002034471357591666031391925350591524680874452002139016807558e-155), SC_(1.86458518280005168582274132886573345934411788365010172356788e-155) }},
            {{ SC_(1.0),  SC_(3.7291703656001033716454826577314669186882357673002034471357591666031391925350591524680874452002139016807558e-155), SC_(1.86458518280005168582274132886573345934411788365010172356788e-155) }},
            {{ SC_(-1.125), SC_(3.7291703656001033716454826577314669186882357673002034471357591666031391925350591524680874452002139016807558e-155), SC_(-1.34963720853101363690381585556234820027343435206156667634081e173) }},
            {{ SC_(1.125),  SC_(3.7291703656001033716454826577314669186882357673002034471357591666031391925350591524680874452002139016807558e-155), SC_(8.02269390325932403421158766283366891170783955777638875887348e-175) }},
            {{ SC_(0.5), SC_(1.2458993688871959419388378518880931736878259938089494331010226962863582408064841833232475731084062642684629e-206), SC_(8.90597649117647254543282704099383321071493400182381039079219e-104) }},
        }};

#include "math/test/bessel_i_int_data.ipp"
#include "math/test/bessel_i_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0], datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
                if (!(actual == datum[2])) {
                    static constexpr std::array<T, 3> bad_value = {{
                        SC_(0.7e2), SC_(0.177219114266335964202880859375e-2), SC_(0.175887342640394106189151976112543057962e-313)
                    }};
                    if (!(::equal(bad_value, datum))) { // TRANSITION, VSO#FIXME
                        BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
                    }
                }
            };
        };

        ::for_each(i0_data, tester(eps<T>));
        ::for_each(i1_data, tester(eps<T>));
        ::for_each(in_data, tester(4 * eps<T>));
        ::for_each(iv_data, tester(4 * eps<T>));
        if (0 != static_cast<T>(std::ldexp(0.5, -700))) {
            ::for_each(iv_large_data, tester(4 * eps<T>));
        }
        ::for_each(bessel_i_int_data, tester(10 * eps<T>));
        ::for_each(bessel_i_data, tester(8 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_i_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace cyl_bessel_i

namespace cyl_bessel_j {
    template<class T>
    constexpr auto control_fn = [](T nu, T x) {
        return boost::math::cyl_bessel_j(nu, x);
    };

    template<class T>
    constexpr auto test_fn = [](auto nu, T) {
        static_assert(always_false<decltype(nu)>);
    };
    template<>
    constexpr auto test_fn<float> = std::cyl_bessel_jf;
    template<>
    constexpr auto test_fn<double> = std::cyl_bessel_j;
    template<>
    constexpr auto test_fn<long double> = std::cyl_bessel_jl;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_j, T, fptypes) {
        // Data taken from test_bessel_j.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 8> j0_data = {{
            {{ SC_(0.0), SC_(0.0), SC_(1.0) }},
            {{ SC_(0.0), SC_(1.0), SC_(0.7651976865579665514497175261026632209093) }},
            {{ SC_(0.0), SC_(-2.0), SC_(0.2238907791412356680518274546499486258252) }},
            {{ SC_(0.0), SC_(4.0), SC_(-0.3971498098638473722865907684516980419756) }},
            {{ SC_(0.0), SC_(-8.0), SC_(0.1716508071375539060908694078519720010684) }},
            {{ SC_(0.0), SC_(1e-05), SC_(0.999999999975000000000156249999999565972) }},
            {{ SC_(0.0), SC_(1e-10), SC_(0.999999999999999999997500000000000000000) }},
            {{ SC_(0.0), SC_(-1e+01), SC_(-0.2459357644513483351977608624853287538296) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 6> j0_tricky = {{
            // Big numbers make the accuracy of std::sin the limiting factor:
            {{ SC_(0.0), SC_(1e+03), SC_(0.02478668615242017456133073111569370878617) }},
            {{ SC_(0.0), SC_(1e+05), SC_(-0.001719201116235972192570601477073201747532) }},
            // test at the roots:
            {{ SC_(0.0), SC_(2.4048252105712890625) /*T(2521642.0) / (1024 * 1024)*/, SC_(1.80208819970046790002973759410972422387259992955354630042138e-7) }},
            {{ SC_(0.0), SC_(5.52007770538330078125) /*T(5788221.0) / (1024 * 1024)*/, SC_(-1.37774249380686777043369399806210229535671843632174587432454e-7) }},
            {{ SC_(0.0), SC_(8.65372753143310546875) /*T(9074091.0) / (1024 * 1024)*/, SC_(1.03553057441100845081018471279571355857520645127532785991335e-7) }},
            {{ SC_(0.0), SC_(11.791534423828125) /*T(12364320.0) / (1024 * 1024)*/, SC_(-3.53017140778223781420794006033810387155048392363051866610931e-9) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 8> j1_data = {{
            {{ SC_(1.0), SC_(0.0), SC_(0.0) }},
            {{ SC_(1.0), SC_(1.0), SC_(0.4400505857449335159596822037189149131274) }},
            {{ SC_(1.0), SC_(-2.0), SC_(-0.5767248077568733872024482422691370869203) }},
            {{ SC_(1.0), SC_(4.0), SC_(-6.604332802354913614318542080327502872742e-02) }},
            {{ SC_(1.0), SC_(-8.0), SC_(-0.2346363468539146243812766515904546115488) }},
            {{ SC_(1.0), SC_(1e-05), SC_(4.999999999937500000000260416666666124132e-06) }},
            {{ SC_(1.0), SC_(1e-10), SC_(4.999999999999999999993750000000000000000e-11) }},
            {{ SC_(1.0), SC_(-1e+01), SC_(-4.347274616886143666974876802585928830627e-02) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 5> j1_tricky = {{
            // Big numbers make the accuracy of std::sin the limiting factor:
            {{ SC_(1.0), SC_(1e+03), SC_(4.728311907089523917576071901216916285418e-03) }},
            {{ SC_(1.0), SC_(1e+05), SC_(1.846757562882567716362123967114215743694e-03) }},
            // test zeros:
            {{ SC_(1.0), SC_(3.8317050933837890625) /*T(4017834) / (1024 * 1024)*/, SC_(3.53149033321258645807835062770856949751958513973522222203044e-7) }},
            {{ SC_(1.0), SC_(7.01558589935302734375) /*T(7356375) / (1024 * 1024)*/, SC_(-2.31227973111067286051984021150135526024117175836722748404342e-7) }},
            {{ SC_(1.0), SC_(10.1734676361083984375) /*T(10667654) / (1024 * 1024)*/, SC_(1.24591331097191900488116495350277530373473085499043086981229e-7) }},
        }};

        static const boost::array<boost::array<typename table_type<T>::type, 3>, 17> jn_data = {{
            // This first one is a modified test case from https://svn.boost.org/trac/boost/ticket/2733
            {{ SC_(-1.0), SC_(1.25), SC_(-0.510623260319880467069474837274910375352924050139633057168856) }},
            {{ SC_(2.0), SC_(0.0), SC_(0.0) }},
            {{ SC_(-2.0), SC_(0.0), SC_(0.0) }},
            {{ SC_(2.0), SC_(1e-02), SC_(1.249989583365885362413250958437642113452e-05) }},
            {{ SC_(5.0), SC_(10.0), SC_(-0.2340615281867936404436949416457777864635) }},
            {{ SC_(5.0), SC_(-10.0), SC_(0.2340615281867936404436949416457777864635) }},
            {{ SC_(-5.0), SC_(1e+06), SC_(7.259643842453285052375779970433848914846e-04) }},
            {{ SC_(5.0), SC_(1e+06), SC_(-0.000725964384245328505237577997043384891484649290328285235308619) }},
            {{ SC_(-5.0), SC_(-1.0), SC_(2.497577302112344313750655409880451981584e-04) }},
            {{ SC_(10.0), SC_(10.0), SC_(0.2074861066333588576972787235187534280327) }},
            {{ SC_(10.0), SC_(-10.0), SC_(0.2074861066333588576972787235187534280327) }},
            {{ SC_(10.0), SC_(-5.0), SC_(1.467802647310474131107532232606627020895e-03) }},
            {{ SC_(-10.0), SC_(1e+06), SC_(-3.310793117604488741264958559035744460210e-04) }},
            {{ SC_(10.0), SC_(1e+06), SC_(-0.000331079311760448874126495855903574446020957243277028930713243) }},
            {{ SC_(1e+02), SC_(8e+01), SC_(4.606553064823477354141298259169874909670e-06) }},
            {{ SC_(1e+03), SC_(1e+05), SC_(1.283178112502480365195139312635384057363e-03) }},
            {{ SC_(10.0), SC_(1e-100), SC_(2.69114445546737213403880070546737213403880070546737213403880e-1010) }},
        }};

        static const boost::array<boost::array<typename table_type<T>::type, 3>, 20> jv_data = {{
            //SC_(-2.4), {{ SC_(0.0), std::numeric_limits<T>::infinity() }},
            {{ SC_(22.5), SC_(0.0), SC_(0.0) }},
            {{ SC_(2.3994140625) /*2457.0 / 1024*/, SC_(0.0009765625) /* 1 / 1024*/, SC_(3.80739920118603335646474073457326714709615200130620574875292e-9) }},
            {{ SC_(5.5), SC_(3.1416015625) /* 3217/1024*/, SC_(0.0281933076257506091621579544064767140470089107926550720453038) }},
            {{ SC_(-5.5), SC_(3.1416015625) /* 3217/1024*/, SC_(-2.55820064470647911823175836997490971806135336759164272675969) }},
            {{ SC_(-5.5), SC_(1e+04), SC_(2.449843111985605522111159013846599118397e-03) }},
            {{ SC_(5.5), SC_(1e+04), SC_(0.00759343502722670361395585198154817047185480147294665270646578) }},
            {{ SC_(5.5), SC_(1e+06), SC_(-0.000747424248595630177396350688505919533097973148718960064663632) }},
            {{ SC_(5.125), SC_(1e+06), SC_(-0.000776600124835704280633640911329691642748783663198207360238214) }},
            {{ SC_(5.875), SC_(1e+06), SC_(-0.000466322721115193071631008581529503095819705088484386434589780) }},
            {{ SC_(0.5), SC_(101.0), SC_(0.0358874487875643822020496677692429287863419555699447066226409) }},
            {{ SC_(-5.5), SC_(1e+04), SC_(0.00244984311198560552211115901384659911839737686676766460822577) }},
            {{ SC_(-5.5), SC_(1e+06), SC_(0.000279243200433579511095229508894156656558211060453622750659554) }},
            {{ SC_(-0.5), SC_(101.0), SC_(0.0708184798097594268482290389188138201440114881159344944791454) }},
            {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(0.0009765625) /* 1/1024*/, SC_(1.41474013160494695750009004222225969090304185981836460288562e35) }},
            {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(15.0), SC_(-0.0902239288885423309568944543848111461724911781719692852541489) }},
            {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(100.0), SC_(-0.05476136603168065513386371539426045507795139476742228638) }},
            {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(20000.0), SC_(-0.00556869085445857782456414284057389040183758546505700058) }},
            // Bug report https://svn.boost.org/trac/boost/ticket/4812:
            {{ SC_(1.5), SC_(7.845703125) /* 8034/1024*/, SC_(0.0339477646369710610146236955872928005087352629422508823945264) }},
            {{ SC_(8.5), SC_(12.566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271928592346053) /*Pi * 4*/, SC_(0.0436807946352780974532519564114026730332781693877984686758680) }},
            {{ SC_(-8.5), SC_(12.566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271928592346053) /*Pi * 4*/, SC_(-0.257086543428224355151772807588810984369026142375675714560864) }},
        }};

        static const boost::array<boost::array<typename table_type<T>::type, 3>, 4> jv_large_data = {{
            // Bug report https://svn.boost.org/trac/boost/ticket/5560:
            {{ SC_(-0.5), SC_(1.2458993688871959419388378518880931736878259938089494331010226962863582408064841833232475731084062642684629e-206) /*static_cast<T>(std::ldexp(0.5, -683))*/, SC_(7.14823099969225685526188875418476476336424046896822867989728e102) }},
            {{ SC_(256.0), SC_(512.0), SC_(0.00671672065717513246956991122723250578101154313313749938944675) }},
            {{ SC_(-256.0), SC_(8.0), SC_(1.46866142030022704638298523775638527553596432641223316232692e-353) }},
            {{ SC_(-2.5), SC_(4.0), SC_(-0.0145679476685218007666785535204236327832335803441449596297004) }},
        }};


#include "math/test/bessel_j_int_data.ipp"
#include "math/test/bessel_j_data.ipp"
#include "math/test/bessel_j_large_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0], datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
                if (!(actual == datum[2])) {
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
                }
            };
        };

        ::for_each(j0_data, tester(4 * eps<T>));
        ::for_each(j0_tricky, tester(10030000 * eps<T>));
        ::for_each(j1_data, tester(12 * eps<T>));
        ::for_each(j1_tricky, tester(35000 * eps<T>));
        ::for_each(jn_data, tester(15 * eps<T>));
        ::for_each(jv_data, tester(20 * eps<T>));
        if (static_cast<T>(jv_large_data[0][1]) != 0)
            ::for_each(jv_large_data, tester(10 * eps<T>));
        ::for_each(bessel_j_int_data, tester(20 * eps<T>));
        ::for_each(bessel_j_data, tester(10 * eps<T>));
        ::for_each(bessel_j_large_data, tester(60 * eps<T>));

        //
        // Some special cases:
        //
        BOOST_CHECK_EQUAL(test_fn<T>(0, T(0)), T(1));
        BOOST_CHECK_EQUAL(test_fn<T>(1, T(0)), T(0));
        BOOST_CHECK_EQUAL(test_fn<T>(100000, T(0)), T(0));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_j_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is x >= 0

        //
        // Special cases that are errors:
        //
        BOOST_CHECK(std::isnan(test_fn<T>(T(-2.5), T(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(-2.5), T(-2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(2.5), T(-2))));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace cyl_bessel_j

namespace cyl_bessel_k {
    template<class T>
    constexpr auto control_fn = [](T nu, T x) {
        return boost::math::cyl_bessel_k(nu, x);
    };

    template<class T>
    constexpr auto test_fn = [](auto nu, T) {
        static_assert(always_false<decltype(nu)>);
    };
    template<>
    constexpr auto test_fn<float> = std::cyl_bessel_kf;
    template<>
    constexpr auto test_fn<double> = std::cyl_bessel_k;
    template<>
    constexpr auto test_fn<long double> = std::cyl_bessel_kl;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_k, T, fptypes) {
        // Data taken from test_bessel_k.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 9> k0_data = {{
            {{ SC_(0.0), SC_(1.0), SC_(0.421024438240708333335627379212609036136219748226660472298970) }},
            {{ SC_(0.0), SC_(2.0), SC_(0.113893872749533435652719574932481832998326624388808882892530) }},
            {{ SC_(0.0), SC_(4.0), SC_(0.0111596760858530242697451959798334892250090238884743405382553) }},
            {{ SC_(0.0), SC_(8.0), SC_(0.000146470705222815387096584408698677921967305368833759024089154) }},
            {{ SC_(0.0), SC_(0.000030517578125) /*T(std::ldexp(1.0, -15))*/, SC_(10.5131392267382037062459525561594822400447325776672021972753) }},
            {{ SC_(0.0), SC_(9.31322574615478515625e-10) /*T(std::ldexp(1.0, -30))*/, SC_(20.9103469324567717360787328239372191382743831365906131108531) }},
            {{ SC_(0.0), SC_(8.67361737988403547205962240695953369140625e-19) /*T(std::ldexp(1.0, -60))*/, SC_(41.7047623492551310138446473188663682295952219631968830346918) }},
            {{ SC_(0.0), SC_(50.0), SC_(3.41016774978949551392067551235295223184502537762334808993276e-23) }},
            {{ SC_(0.0), SC_(100.0), SC_(4.65662822917590201893900528948388635580753948544211387402671e-45) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 9> k1_data = {{
            {{ SC_(1.0), SC_(1.0), SC_(0.601907230197234574737540001535617339261586889968106456017768) }},
            {{ SC_(1.0), SC_(2.0), SC_(0.139865881816522427284598807035411023887234584841515530384442) }},
            {{ SC_(1.0), SC_(4.0), SC_(0.0124834988872684314703841799808060684838415849886258457917076) }},
            {{ SC_(1.0), SC_(8.0), SC_(0.000155369211805001133916862450622474621117065122872616157079566) }},
            {{ SC_(1.0), SC_(0.000030517578125) /*T(std::ldexp(1.0, -15))*/, SC_(32767.9998319528316432647441316539139725104728341577594326513) }},
            {{ SC_(1.0), SC_(9.31322574615478515625e-10) /*T(std::ldexp(1.0, -30))*/, SC_(1.07374182399999999003003028572687332810353799544215073362305e9) }},
            {{ SC_(1.0), SC_(8.67361737988403547205962240695953369140625e-19) /*T(std::ldexp(1.0, -60))*/, SC_(1.15292150460684697599999999999999998169660198868126604634036e18) }},
            {{ SC_(1.0), SC_(50.0), SC_(3.44410222671755561259185303591267155099677251348256880221927e-23) }},
            {{ SC_(1.0), SC_(100.0), SC_(4.67985373563690928656254424202433530797494354694335352937465e-45) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 9> kn_data = {{
            {{ SC_(2.0), SC_(9.31322574615478515625e-10) /*T(std::ldexp(1.0, -30))*/, SC_(2.30584300921369395150000000000000000234841952009593636868109e18) }},
            {{ SC_(5.0), SC_(10.0), SC_(0.0000575418499853122792763740236992723196597629124356739596921536) }},
            {{ SC_(-5.0), SC_(100.0), SC_(5.27325611329294989461777188449044716451716555009882448801072e-45) }},
            {{ SC_(10.0), SC_(10.0), SC_(0.00161425530039067002345725193091329085443750382929208307802221) }},
            {{ SC_(10.0), SC_(9.31322574615478515625e-10) /*T(std::ldexp(1.0, -30))*/, SC_(3.78470202927236255215249281534478864916684072926050665209083e98) }},
            {{ SC_(-10.0), SC_(1.0), SC_(1.80713289901029454691597861302340015908245782948536080022119e8) }},
            {{ SC_(100.0), SC_(5.0), SC_(7.03986019306167654653386616796116726248616158936088056952477e115) }},
            {{ SC_(100.0), SC_(80.0), SC_(8.39287107246490782848985384895907681748152272748337807033319e-12) }},
            {{ SC_(-1000.0), SC_(700.0), SC_(6.51561979144735818903553852606383312984409361984128221539405e-31) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 11> kv_data = {{
            {{ SC_(0.5), SC_(0.875), SC_(0.558532231646608646115729767013630967055657943463362504577189) }},
            {{ SC_(0.5), SC_(1.125), SC_(0.383621010650189547146769320487006220295290256657827220786527) }},
            {{ SC_(2.25), SC_(9.31322574615478515625e-10) /*T(std::ldexp(1.0, -30))*/, SC_(5.62397392719283271332307799146649700147907612095185712015604e20) }},
            {{ SC_(5.5), SC_(3.1416015625) /*3217/1024*/, SC_(1.30623288775012596319554857587765179889689223531159532808379) }},
            {{ SC_(-5.5), SC_(10.0), SC_(0.0000733045300798502164644836879577484533096239574909573072142667) }},
            {{ SC_(-5.5), SC_(100.0), SC_(5.41274555306792267322084448693957747924412508020839543293369e-45) }},
            {{ SC_(10.0), SC_(0.0009765625) /*1/1024*/, SC_(2.35522579263922076203415803966825431039900000000993410734978e38) }},
            {{ SC_(10.0), SC_(10.0), SC_(0.00161425530039067002345725193091329085443750382929208307802221) }},
            {{ SC_(141.3994140625) /*T(144793)/1024)*/, SC_(100.0), SC_(1.39565245860302528069481472855619216759142225046370312329416e-6) }},
            {{ SC_(141.3994140625) /*T(144793)/1024)*/, SC_(200.0), SC_(9.11950412043225432171915100042647230802198254567007382956336e-68) }},
            {{ SC_(-141.3994140625) /*T(-144793)/1024*/, SC_(50.0), SC_(1.30185229717525025165362673848737761549946548375142378172956e42) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 5> kv_large_data = {{
            // Bug report https://svn.boost.org/trac/boost/ticket/5560:
            {{ SC_(-1.0), SC_(3.729170365600103371645482657731466918688235767300203447135759166603139192535059152468087445200213901e-155) /*static_cast<T>(ldexp(0.5, -512))*/, SC_(2.68156158598851941991480499964116922549587316411847867554471e154) }},
            {{ SC_(1.0),  SC_(3.729170365600103371645482657731466918688235767300203447135759166603139192535059152468087445200213901e-155) /*static_cast<T>(ldexp(0.5, -512))*/, SC_(2.68156158598851941991480499964116922549587316411847867554471e154) }},
            {{ SC_(-1.125), SC_(3.729170365600103371645482657731466918688235767300203447135759166603139192535059152468087445200213901e-155) /*static_cast<T>(ldexp(0.5, -512))*/, SC_(5.53984048006472105611199242328122729730752165907526178753978e173) }},
            {{ SC_(1.125),  SC_(3.729170365600103371645482657731466918688235767300203447135759166603139192535059152468087445200213901e-155) /*static_cast<T>(ldexp(0.5, -512))*/, SC_(5.53984048006472105611199242328122729730752165907526178753978e173) }},
            {{ SC_(0.5),  SC_(1.2458993688871959419388378518880931736878259938089494331010226962863582408064841833232475731084062642684629e-206) /*static_cast<T>(ldexp(0.5, -683))*/, SC_(1.12284149973980088540335945247019177715948513804063794284101e103) }},
        }};

#include "math/test/bessel_k_int_data.ipp"
#include "math/test/bessel_k_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0], datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
                if (!(actual == datum[2])) {
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
                }
            };
        };

        ::for_each(k0_data, tester(eps<T>));
        ::for_each(k1_data, tester(eps<T>));
        ::for_each(kn_data, tester(4 * eps<T>));
        ::for_each(kv_data, tester(5 * eps<T>));
        if (0 != static_cast<T>(std::ldexp(0.5, -512))) {
            ::for_each(kv_large_data, tester(60 * eps<T>));
        }
        ::for_each(bessel_k_int_data, tester(10 * eps<T>));
        ::for_each(bessel_k_data, tester(10 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_k_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is x >= 0

        //
        // Special cases that are errors:
        //
        BOOST_CHECK(std::isnan(test_fn<T>(T(-2.5), T(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(-2.5), T(-2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(2.5), T(-2))));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace cyl_bessel_k

namespace cyl_neumann {
    template<class T>
    constexpr auto control_fn = [](T nu, T x) {
        return boost::math::cyl_neumann(nu, x);
    };

    template<class T>
    constexpr auto test_fn = [](auto nu, T) {
        static_assert(always_false<decltype(nu)>);
    };
    template<>
    constexpr auto test_fn<float> = std::cyl_neumannf;
    template<>
    constexpr auto test_fn<double> = std::cyl_neumann;
    template<>
    constexpr auto test_fn<long double> = std::cyl_neumannl;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_neumann, T, fptypes) {
        // Data taken from test_bessel_y.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 9> y0_data = {{
            {{ SC_(0.0), SC_(1.0), SC_(0.0882569642156769579829267660235151628278175230906755467110438) }},
            {{ SC_(0.0), SC_(2.0), SC_(0.510375672649745119596606592727157873268139227085846135571839) }},
            {{ SC_(0.0), SC_(4.0), SC_(-0.0169407393250649919036351344471532182404925898980149027169321) }},
            {{ SC_(0.0), SC_(8.0), SC_(0.223521489387566220527323400498620359274814930781423577578334) }},
            {{ SC_(0.0), SC_(1e-05), SC_(-7.40316028370197013259676050746759072070960287586102867247159) }},
            {{ SC_(0.0), SC_(1e-10), SC_(-14.7325162726972420426916696426209144888762342592762415255386) }},
            {{ SC_(0.0), SC_(1e-20), SC_(-29.3912282502857968601858410375186700783698345615477536431464) }},
            {{ SC_(0.0), SC_(1e+03), SC_(0.00471591797762281339977326146566525500985900489680197718528000) }},
            {{ SC_(0.0), SC_(1e+05), SC_(0.00184676615886506410434074102431546125884886798090392516843524) }}
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 9> y1_data = {{
            {{ SC_(1.0), SC_(1.0), SC_(-0.781212821300288716547150000047964820549906390716444607843833) }},
            {{ SC_(1.0), SC_(2.0), SC_(-0.107032431540937546888370772277476636687480898235053860525795) }},
            {{ SC_(1.0), SC_(4.0), SC_(0.397925710557100005253979972450791852271189181622908340876586) }},
            {{ SC_(1.0), SC_(8.0), SC_(-0.158060461731247494255555266187483550327344049526705737651263) }},
            {{ SC_(1.0), SC_(1e-10), SC_(-6.36619772367581343150789184284462611709080831190542841855708e9) }},
            {{ SC_(1.0), SC_(1e-20), SC_(-6.36619772367581343075535053490057448139324059868649274367256e19) }},
            {{ SC_(1.0), SC_(1e+01), SC_(0.249015424206953883923283474663222803260416543069658461246944) }},
            {{ SC_(1.0), SC_(1e+03), SC_(-0.0247843312923517789148623560971412909386318548648705287583490) }},
            {{ SC_(1.0), SC_(1e+05), SC_(0.00171921035008825630099494523539897102954509504993494957572726) }}
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 10> yn_data = {{
            {{ SC_(2.0), SC_(1e-20), SC_(-1.27323954473516268615107010698011489627570899691226996904849e40) }},
            {{ SC_(5.0), SC_(10.0), SC_(0.135403047689362303197029014762241709088405766746419538495983) }},
            {{ SC_(-5.0), SC_(1e+06), SC_(0.000331052088322609048503535570014688967096938338061796192422114) }},
            {{ SC_(10.0), SC_(10.0), SC_(-0.359814152183402722051986577343560609358382147846904467526222) }},
            {{ SC_(10.0), SC_(1e-10), SC_(-1.18280490494334933900960937719565669877576135140014365217993e108) }},
            {{ SC_(-10.0), SC_(1e+06), SC_(0.000725951969295187086245251366365393653610914686201194434805730) }},
            {{ SC_(1e+02), SC_(5.0), SC_(-5.08486391602022287993091563093082035595081274976837280338134e115) }},
            {{ SC_(1e+03), SC_(1e+05), SC_(0.00217254919137684037092834146629212647764581965821326561261181) }},
            {{ SC_(-1e+03), SC_(7e+02), SC_(-1.88753109980945889960843803284345261796244752396992106755091e77) }},
            {{ SC_(-25.0), SC_(8.0), SC_(3.45113613777297661997458045843868931827873456761831907587263e8) }}
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 11> yv_data = {{
            //SC_(2.25), {{ SC_(1.0) / 1024, SC_(-1.01759203636941035147948317764932151601257765988969544340275e7) }},
            {{ SC_(0.5), SC_(9.5367431640625e-7) /* 1/(1024*1024)*/, SC_(-817.033790261762580469303126467917092806755460418223776544122) }},
            {{ SC_(5.5), SC_(3.125), SC_(-2.61489440328417468776474188539366752698192046890955453259866) }},
            {{ SC_(-5.5), SC_(3.125), SC_(-0.0274994493896489729948109971802244976377957234563871795364056) }},
            {{ SC_(-5.5), SC_(1e+04), SC_(-0.00759343502722670361395585198154817047185480147294665270646578) }},
            {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(0.0009765625) /*1/1024*/, SC_(-1.50382374389531766117868938966858995093408410498915220070230e38) }},
            {{ SC_(-10.0002994537353515625) /* -10486074 / (1024*1024)*/, SC_(1e+02), SC_(0.0583041891319026009955779707640455341990844522293730214223545) }},
            {{ SC_(141.75), SC_(1e+02), SC_(-5.38829231428696507293191118661269920130838607482708483122068e9) }},
            {{ SC_(141.75), SC_(2e+04), SC_(-0.00376577888677186194728129112270988602876597726657372330194186) }},
            {{ SC_(-141.75), SC_(1e+02), SC_(-3.81009803444766877495905954105669819951653361036342457919021e9) }},
            {{ SC_(8.5), SC_(12.56637061435917295385057353311801153678867759750042328389) /*4Pi*/, SC_(0.257086543428224355151772807588810984369026142375675714560864) }},
            {{ SC_(-8.5), SC_(12.56637061435917295385057353311801153678867759750042328389) /*4Pi*/, SC_(0.0436807946352780974532519564114026730332781693877984686758680) }},
        }};
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 7> yv_large_data = {{
            // Bug report https://svn.boost.org/trac/boost/ticket/5560:
            {{ SC_(0.5), SC_(1.24589936888719594193883785188809317368782599380894e-206) /*static_cast<T>(std::ldexp(0.5, -683))*/, SC_(-7.14823099969225685526188875418476476336424046896822867989728e102) }},
            {{ SC_(-0.5), SC_(1.24589936888719594193883785188809317368782599380894e-206) /*static_cast<T>(std::ldexp(0.5, -683))*/, SC_(8.90597649117647254543282704099383321071493400182381039079219e-104) }},
            {{ SC_(0.0), SC_(1.1102230246251565404236316680908203125e-16) /*static_cast<T>(std::ldexp(1.0, -53))*/, SC_(-23.4611779112897561252987257324561640034037313549011724328997) }},
            {{ SC_(1.0), SC_(1.1102230246251565404236316680908203125e-16) /*static_cast<T>(std::ldexp(1.0, -53))*/, SC_(-5.73416113922265864550047623401604244038331542638719289100990e15) }},
            {{ SC_(2.0), SC_(1.1102230246251565404236316680908203125e-16) /*static_cast<T>(std::ldexp(1.0, -53))*/, SC_(-1.03297463879542177245046832533417970379386617249046560049244e32) }},
            {{ SC_(3.0), SC_(1.1102230246251565404236316680908203125e-16) /*static_cast<T>(std::ldexp(1.0, -53))*/, SC_(-3.72168335868978735639260528876490232745489151562358712422544e48) }},
            {{ SC_(10.0), SC_(1.1102230246251565404236316680908203125e-16) /*static_cast<T>(std::ldexp(1.0, -53))*/, SC_(-4.15729476804920974669173904282420477878640623992500096231384e167) }},
        }};

#include "math/test/bessel_y01_data.ipp"
#include "math/test/bessel_yn_data.ipp"
#include "math/test/bessel_yv_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0], datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0], datum[1]));
                if (!(actual == datum[2])) {
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
                }
            };
        };

        ::for_each(y0_data, tester(5 * eps<T>));
        ::for_each(y1_data, tester(5 * eps<T>));
        ::for_each(yn_data, tester(35 * eps<T>));
        ::for_each(yv_data, tester(14 * eps<T>));
        if (static_cast<T>(yv_large_data[0][1]) != 0) {
            ::for_each(yv_large_data, tester(eps<T>));
        }
        ::for_each(bessel_y01_data, tester(5 * eps<T>));
        ::for_each(bessel_yn_data, tester(120 * eps<T>));
        ::for_each(bessel_yv_data, tester(1250 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_neumann_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is x >= 0

        //
        // Special cases that are errors:
        //
        BOOST_CHECK(std::isnan(test_fn<T>(T(-2.5), T(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(-2.5), T(-2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(2.5), T(-2))));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace cyl_neumann

namespace ellint_1 {
    template<class T>
    constexpr auto control_fn = [](T k, T phi) {
        return boost::math::ellint_1(k, phi);
    };

    template<class T>
    constexpr auto test_fn = [](auto k, T) {
        static_assert(always_false<decltype(k)>);
    };
    template<>
    constexpr auto test_fn<float> = std::ellint_1f;
    template<>
    constexpr auto test_fn<double> = std::ellint_1;
    template<>
    constexpr auto test_fn<long double> = std::ellint_1l;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_ellint_1, T, fptypes) {
        // From function test_spots in test_ellint_1.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 19> data1 = {{
            {{ SC_(0.0), SC_(0.0), SC_(0.0) }},
            {{ SC_(-10.0), SC_(0.0), SC_(-10.0) }},
            {{ SC_(-1.0), SC_(-1.0), SC_(-1.2261911708835170708130609674719067527242483502207) }},
            {{ SC_(-4.0), SC_(0.875), SC_(-5.3190556182262405182189463092940736859067548232647) }},
            {{ SC_(8.0), SC_(-0.625), SC_(9.0419973860310100524448893214394562615252527557062) }},
            {{ SC_(1e-05), SC_(0.875), SC_(0.000010000000000127604166668510945638036143355898993088) }},
            {{ SC_(1e+05), SC_(0.009765625) /*T(10)/1024*/, SC_(100002.38431454899771096037307519328741455615271038) }},
            {{ SC_(1e-20), SC_(1.0), SC_(1.0000000000000000000000000000000000000000166666667e-20) }},
            {{ SC_(1e-20), SC_(1e-20), SC_(1.000000000000000e-20) }},
            {{ SC_(1e+20), SC_(0.390625) /*T(400)/1024*/, SC_(1.0418143796499216839719289963154558027005142709763e20) }},
            {{ SC_(1e+50), SC_(0.875), SC_(1.3913251718238765549409892714295358043696028445944e50) }},
            {{ SC_(2.0), SC_(0.5), SC_(2.1765877052210673672479877957388515321497888026770) }},
            {{ SC_(4.0), SC_(0.5), SC_(4.2543274975235836861894752787874633017836785640477) }},
            {{ SC_(6.0), SC_(0.5), SC_(6.4588766202317746302999080620490579800463614807916) }},
            {{ SC_(10.0), SC_(0.5), SC_(10.697409951222544858346795279378531495869386960090) }},
            {{ SC_(-2.0), SC_(0.5), SC_(-2.1765877052210673672479877957388515321497888026770) }},
            {{ SC_(-4.0), SC_(0.5), SC_(-4.2543274975235836861894752787874633017836785640477) }},
            {{ SC_(-6.0), SC_(0.5), SC_(-6.4588766202317746302999080620490579800463614807916) }},
            {{ SC_(-10.0), SC_(0.5), SC_(-10.697409951222544858346795279378531495869386960090) }},
        }};

#include "math/test/ellint_f_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[1], datum[0]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[1], datum[0]));
                if (!(actual == datum[2])) // +/-inf is equal to, but not "close" to, +/-inf
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
            };
        };

        ::for_each(data1, tester(eps<T>));
        ::for_each(ellint_f_data, tester(3 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_ellint_1_boundaries, T, fptypes) {
        auto const tolerance = eps<T>;

        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |k| <= 1
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(2), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(-2), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1), static_cast<T>(0)),
            static_cast<T>(0L), tolerance);
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(0), static_cast<T>(1)),
            static_cast<T>(1L), tolerance);
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace ellint_1

namespace ellint_2 {
    template<class T>
    constexpr auto control_fn = [](T k, T phi) {
        return boost::math::ellint_2(k, phi);
    };

    template<class T>
    constexpr auto test_fn = [](auto k, T) {
        static_assert(always_false<decltype(k)>);
    };
    template<>
    constexpr auto test_fn<float> = std::ellint_2f;
    template<>
    constexpr auto test_fn<double> = std::ellint_2;
    template<>
    constexpr auto test_fn<long double> = std::ellint_2l;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_ellint_2, T, fptypes) {
        // From function test_spots in test_ellint_2.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 3>, 10> data1 = {{
            {{ SC_(0.0), SC_(0.0), SC_(0.0) }},
            {{ SC_(-10.0), SC_(0.0), SC_(-10.0) }},
            {{ SC_(-1.0), SC_(-1.0), SC_(-0.84147098480789650665250232163029899962256306079837) }},
            {{ SC_(-4.0), SC_(0.87890625) /*T(900) / 1024*/, SC_(-3.1756145986492562317862928524528520686391383168377) }},
            {{ SC_(8.0), SC_(-0.5859375) /*T(-600) / 1024*/, SC_(7.2473147180505693037677015377802777959345489333465) }},
            {{ SC_(1e-05), SC_(0.78125) /*T(800) / 1024*/, SC_(9.999999999898274739584436515967055859383969942432E-6) }},
            {{ SC_(1e+05), SC_(0.09765625) /*T(100) / 1024*/, SC_(99761.153306972066658135668386691227343323331995888) }},
            {{ SC_(1e+10), SC_(-0.5), SC_(9.3421545766487137036576748555295222252286528414669e9) }},
            {{ SC_(7.3786976294838206464e19) /*static_cast<T>(ldexp(T(1), 66))*/, SC_(0.390625) /*T(400) / 1024*/, SC_(7.0886102721911705466476846969992069994308167515242e19) }},
            {{ SC_(9.3536104789177786765035829293842113257979682750464e49) /*static_cast<T>(ldexp(T(1), 166))*/, SC_(0.87890625) /*T(900) / 1024*/, SC_(7.1259011068364515942912094521783688927118026465790e49) }},
        }};

#include "math/test/ellint_e2_data.ipp"

        static const boost::array<boost::array<typename table_type<T>::type, 3>, 72> small_angles = { {
            {{ SC_(0.00097656250000000000000000000000000000000000000000000), SC_(0.5), SC_(0.00097656246119489873806295171767681128826061680891539) }},{{ SC_(0.00048828125000000000000000000000000000000000000000000), SC_(0.5), SC_(0.00048828124514936177847275804383491089917731651869089) }},{{ SC_(0.00024414062500000000000000000000000000000000000000000), SC_(0.5), SC_(0.00024414062439367020469080959924292294147407037569089) }},{{ SC_(0.00012207031250000000000000000000000000000000000000000), SC_(0.5), SC_(0.00012207031242420877503577978533579671450656676021144) }},{{ SC_(0.000061035156250000000000000000000000000000000000000000), SC_(0.5), SC_(0.000061035156240526096862267116434822602398026203555135) }},{{ SC_(0.000030517578125000000000000000000000000000000000000000), SC_(0.5), SC_(0.000030517578123815762107245722156263286910312978330942) }},{{ SC_(0.000015258789062500000000000000000000000000000000000000), SC_(0.5), SC_(0.000015258789062351970263388913163340973814136929083865) }},{{ SC_(7.6293945312500000000000000000000000000000000000000e-6), SC_(0.5), SC_(7.6293945312314962829230890795991108894734462126936e-6) }},{{ SC_(3.8146972656250000000000000000000000000000000000000e-6), SC_(0.5), SC_(3.8146972656226870353653697266430602974836595561218e-6) }},{{ SC_(1.9073486328125000000000000000000000000000000000000e-6), SC_(0.5), SC_(1.9073486328122108794206707030707941437882919842603e-6) }},{{ SC_(9.5367431640625000000000000000000000000000000000000e-7), SC_(0.5), SC_(9.5367431640621385992758382186011213067376941980499e-7) }},{{ SC_(4.7683715820312500000000000000000000000000000000000e-7), SC_(0.5), SC_(4.7683715820312048249094797723177223079355575583105e-7) }},{{ SC_(2.3841857910156250000000000000000000000000000000000e-7), SC_(0.5), SC_(2.3841857910156193531136849713832334805104830239272e-7) }},{{ SC_(1.1920928955078125000000000000000000000000000000000e-7), SC_(0.5), SC_(1.1920928955078117941392106214180141285643896716624e-7) }},{{ SC_(5.9604644775390625000000000000000000000000000000000e-8), SC_(0.5), SC_(5.9604644775390616176740132767709895180494181165759e-8) }},{{ SC_(2.9802322387695312500000000000000000000000000000000e-8), SC_(0.5), SC_(2.9802322387695311397092516595963259352981751091479e-8) }},{{ SC_(1.4901161193847656250000000000000000000000000000000e-8), SC_(0.5), SC_(1.4901161193847656112136564574495392495854593212863e-8) }},{{ SC_(7.4505805969238281250000000000000000000000000000000e-9), SC_(0.5), SC_(7.4505805969238281077670705718119235956296952243088e-9) }},{{ SC_(3.7252902984619140625000000000000000000000000000000e-9), SC_(0.5), SC_(3.7252902984619140603458838214764904348802078740605e-9) }},{{ SC_(1.8626451492309570312500000000000000000000000000000e-9), SC_(0.5), SC_(1.8626451492309570309807354776845613039046039833520e-9) }},{{ SC_(9.3132257461547851562500000000000000000000000000000e-10), SC_(0.5), SC_(9.3132257461547851559134193471057016297384356039070e-10) }},{{ SC_(4.6566128730773925781250000000000000000000000000000e-10), SC_(0.5), SC_(4.6566128730773925780829274183882127037128569700108e-10) }},{{ SC_(2.3283064365386962890625000000000000000000000000000e-10), SC_(0.5), SC_(2.3283064365386962890572409272985265879639681374864e-10) }},{{ SC_(1.1641532182693481445312500000000000000000000000000e-10), SC_(0.5), SC_(1.1641532182693481445305926159123158234954916739431e-10) }},{{ SC_(5.8207660913467407226562500000000000000000000000000e-11), SC_(0.5), SC_(5.8207660913467407226554282698903947793693632351656e-11) }},{{ SC_(2.9103830456733703613281250000000000000000000000000e-11), SC_(0.5), SC_(2.9103830456733703613280222837362993474211703619812e-11) }},{{ SC_(1.4551915228366851806640625000000000000000000000000e-11), SC_(0.5), SC_(1.4551915228366851806640496604670374184276462939222e-11) }},{{ SC_(7.2759576141834259033203125000000000000000000000000e-12), SC_(0.5), SC_(7.2759576141834259033202964505837967730345578669885e-12) }},{{ SC_(3.6379788070917129516601562500000000000000000000000e-12), SC_(0.5), SC_(3.6379788070917129516601542438229745966293197333606e-12) }},{{ SC_(1.8189894035458564758300781250000000000000000000000e-12), SC_(0.5), SC_(1.8189894035458564758300778742278718245786649666697e-12) }},{{ SC_(9.0949470177292823791503906250000000000000000000000e-13), SC_(0.5), SC_(9.0949470177292823791503903115348397807233312083370e-13) }},{{ SC_(4.5474735088646411895751953125000000000000000000000e-13), SC_(0.5), SC_(4.5474735088646411895751952733168549725904164010421e-13) }},{{ SC_(2.2737367544323205947875976562500000000000000000000e-13), SC_(0.5), SC_(2.2737367544323205947875976513521068715738020501303e-13) }},{{ SC_(1.1368683772161602973937988281250000000000000000000e-13), SC_(0.5), SC_(1.1368683772161602973937988275127633589467252562663e-13) }},{{ SC_(5.6843418860808014869689941406250000000000000000000e-14), SC_(0.5), SC_(5.6843418860808014869689941398597041986834065703329e-14) }},{{ SC_(2.8421709430404007434844970703125000000000000000000e-14), SC_(0.5), SC_(2.8421709430404007434844970702168380248354258212916e-14) }},{{ SC_(1.4210854715202003717422485351562500000000000000000e-14), SC_(0.5), SC_(1.4210854715202003717422485351442922531044282276615e-14) }},{{ SC_(7.1054273576010018587112426757812500000000000000000e-15), SC_(0.5), SC_(7.1054273576010018587112426757663028163805352845768e-15) }},{{ SC_(3.5527136788005009293556213378906250000000000000000e-15), SC_(0.5), SC_(3.5527136788005009293556213378887566020475669105721e-15) }},{{ SC_(1.7763568394002504646778106689453125000000000000000e-15), SC_(0.5), SC_(1.7763568394002504646778106689450789502559458638215e-15) }},{{ SC_(8.8817841970012523233890533447265625000000000000000e-16), SC_(0.5), SC_(8.8817841970012523233890533447262705628199323297769e-16) }},{{ SC_(4.4408920985006261616945266723632812500000000000000e-16), SC_(0.5), SC_(4.4408920985006261616945266723632447578524915412221e-16) }},{{ SC_(2.2204460492503130808472633361816406250000000000000e-16), SC_(0.5), SC_(2.2204460492503130808472633361816360634815614426528e-16) }},{{ SC_(1.1102230246251565404236316680908203125000000000000e-16), SC_(0.5), SC_(1.1102230246251565404236316680908197423101951803316e-16) }},{{ SC_(5.5511151231257827021181583404541015625000000000000e-17), SC_(0.5), SC_(5.5511151231257827021181583404541008497627439754145e-17) }},{{ SC_(2.7755575615628913510590791702270507812500000000000e-17), SC_(0.5), SC_(2.7755575615628913510590791702270506921578429969268e-17) }},{{ SC_(1.3877787807814456755295395851135253906250000000000e-17), SC_(0.5), SC_(1.3877787807814456755295395851135253794884803746159e-17) }},{{ SC_(6.9388939039072283776476979255676269531250000000000e-18), SC_(0.5), SC_(6.9388939039072283776476979255676269392043504682698e-18) }},{{ SC_(3.4694469519536141888238489627838134765625000000000e-18), SC_(0.5), SC_(3.4694469519536141888238489627838134748224188085337e-18) }},{{ SC_(1.7347234759768070944119244813919067382812500000000e-18), SC_(0.5), SC_(1.7347234759768070944119244813919067380637398510667e-18) }},{{ SC_(8.6736173798840354720596224069595336914062500000000e-19), SC_(0.5), SC_(8.6736173798840354720596224069595336911343623138334e-19) }},{{ SC_(4.3368086899420177360298112034797668457031250000000e-19), SC_(0.5), SC_(4.3368086899420177360298112034797668456691390392292e-19) }},{{ SC_(2.1684043449710088680149056017398834228515625000000e-19), SC_(0.5), SC_(2.1684043449710088680149056017398834228473142549036e-19) }},{{ SC_(1.0842021724855044340074528008699417114257812500000e-19), SC_(0.5), SC_(1.0842021724855044340074528008699417114252502193630e-19) }},{{ SC_(5.4210108624275221700372640043497085571289062500000e-20), SC_(0.5), SC_(5.4210108624275221700372640043497085571282424617037e-20) }},{{ SC_(2.7105054312137610850186320021748542785644531250000e-20), SC_(0.5), SC_(2.7105054312137610850186320021748542785643701514630e-20) }},{{ SC_(1.3552527156068805425093160010874271392822265625000e-20), SC_(0.5), SC_(1.3552527156068805425093160010874271392822161908079e-20) }},{{ SC_(6.7762635780344027125465800054371356964111328125000e-21), SC_(0.5), SC_(6.7762635780344027125465800054371356964111198478848e-21) }},{{ SC_(3.3881317890172013562732900027185678482055664062500e-21), SC_(0.5), SC_(3.3881317890172013562732900027185678482055647856731e-21) }},{{ SC_(1.6940658945086006781366450013592839241027832031250e-21), SC_(0.5), SC_(1.6940658945086006781366450013592839241027830005529e-21) }},{{ SC_(8.4703294725430033906832250067964196205139160156250e-22), SC_(0.5), SC_(8.4703294725430033906832250067964196205139157624099e-22) }},{{ SC_(4.2351647362715016953416125033982098102569580078125e-22), SC_(0.5), SC_(4.2351647362715016953416125033982098102569579761606e-22) }},{{ SC_(2.1175823681357508476708062516991049051284790039062e-22), SC_(0.5), SC_(2.1175823681357508476708062516991049051284789999498e-22) }},{{ SC_(1.0587911840678754238354031258495524525642395019531e-22), SC_(0.5), SC_(1.0587911840678754238354031258495524525642395014586e-22) }},{{ SC_(5.2939559203393771191770156292477622628211975097656e-23), SC_(0.5), SC_(5.2939559203393771191770156292477622628211975091474e-23) }},{{ SC_(2.6469779601696885595885078146238811314105987548828e-23), SC_(0.5), SC_(2.6469779601696885595885078146238811314105987548055e-23) }},{{ SC_(1.3234889800848442797942539073119405657052993774414e-23), SC_(0.5), SC_(1.3234889800848442797942539073119405657052993774317e-23) }},{{ SC_(6.6174449004242213989712695365597028285264968872070e-24), SC_(0.5), SC_(6.6174449004242213989712695365597028285264968871950e-24) }},{{ SC_(3.3087224502121106994856347682798514142632484436035e-24), SC_(0.5), SC_(3.3087224502121106994856347682798514142632484436020e-24) }},{{ SC_(1.6543612251060553497428173841399257071316242218018e-24), SC_(0.5), SC_(1.6543612251060553497428173841399257071316242218016e-24) }},{{ SC_(8.2718061255302767487140869206996285356581211090088e-25), SC_(0.5), SC_(8.2718061255302767487140869206996285356581211090086e-25) }},{{ SC_(4.1359030627651383743570434603498142678290605545044e-25), SC_(0.5), SC_(4.1359030627651383743570434603498142678290605545044e-25) }}
        } };

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[1], datum[0]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[1], datum[0]));
                if (!(actual == datum[2])) // +/-inf is equal to, but not "close" to, +/-inf
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
            };
        };

        ::for_each(data1, tester(1.5 * eps<T>));
        ::for_each(ellint_e2_data, tester(2.5 * eps<T>));
        ::for_each(small_angles, tester(1.25 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_ellint_2_boundaries, T, fptypes) {
        auto const tolerance = eps<T>;

        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |k| <= 1
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(2), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(-2), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1), static_cast<T>(0)),
            static_cast<T>(0L), tolerance);
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(0), static_cast<T>(1)),
            static_cast<T>(1L), tolerance);
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace ellint_2

namespace ellint_3 {
    template<class T>
    constexpr auto control_fn = [](T k, T nu, T phi) {
        return boost::math::ellint_3(k, nu, phi);
    };

    template<class T>
    constexpr auto test_fn = [](auto k, T, T) {
        static_assert(always_false<decltype(k)>);
    };
    template<>
    constexpr auto test_fn<float> = std::ellint_3f;
    template<>
    constexpr auto test_fn<double> = std::ellint_3;
    template<>
    constexpr auto test_fn<long double> = std::ellint_3l;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_ellint_3, T, fptypes) {
        // From function test_spots in test_ellint_3.hpp:
        static const boost::array<boost::array<typename table_type<T>::type, 4>, 65> data1 = {{
            {{ SC_(1.0), SC_(-1.0), SC_(0.0), SC_(-1.557407724654902230506974807458360173087) }},
            {{ SC_(0.0), SC_(-4.0), SC_(0.4), SC_(-4.153623371196831087495427530365430979011) }},
            {{ SC_(0.0), SC_(8.0), SC_(-0.6), SC_(8.935930619078575123490612395578518914416) }},
            {{ SC_(0.0), SC_(0.5), SC_(0.25), SC_(0.501246705365439492445236118603525029757890291780157969500480) }},
            {{ SC_(0.0), SC_(0.5), SC_(0.0), SC_(0.5) }},
            {{ SC_(-2.0), SC_(0.5), SC_(0.0), SC_(0.437501067017546278595664813509803743009132067629603474488486) }},
            {{ SC_(0.25), SC_(0.5), SC_(0.0), SC_(0.510269830229213412212501938035914557628394166585442994564135) }},
            {{ SC_(0.75), SC_(0.5), SC_(0.0), SC_(0.533293253875952645421201146925578536430596894471541312806165) }},
            {{ SC_(0.75), SC_(0.75), SC_(0.0), SC_(0.871827580412760575085768367421866079353646112288567703061975) }},
            {{ SC_(1.0), SC_(0.25), SC_(0.0), SC_(0.255341921221036266504482236490473678204201638800822621740476) }},
            {{ SC_(2.0), SC_(0.25), SC_(0.0), SC_(0.261119051639220165094943572468224137699644963125853641716219) }},
            {{ SC_(0.9990234375) /*T(1023)/1024*/, SC_(1.5), SC_(0.0), SC_(13.2821612239764190363647953338544569682942329604483733197131) }},
            {{ SC_(0.5), SC_(-1.0), SC_(0.5), SC_(-1.228014414316220642611298946293865487807) }},
            {{ SC_(0.5), SC_(1e+10), SC_(0.5), SC_(1.536591003599172091573590441336982730551e+10) }},
            {{ SC_(-1e+05), SC_(10.0), SC_(0.75), SC_(0.0347926099493147087821620459290460547131012904008557007934290) }},
            {{ SC_(-1e+10), SC_(10.0), SC_(0.875), SC_(0.000109956202759561502329123384755016959364346382187364656768212) }},
            {{ SC_(-1e+10), SC_(1e+20), SC_(0.875), SC_(1.00000626665567332602765201107198822183913978895904937646809e15) }},
            {{ SC_(-1e+10), SC_(1.5703125) /*T(1608)/1024*/, SC_(0.875), SC_(0.0000157080616044072676127333183571107873332593142625043567690379) }},
            {{ SC_(0.9990234375) /*1-T(1) / 1024*/, SC_(1e+20), SC_(0.875), SC_(6.43274293944380717581167058274600202023334985100499739678963e21) }},
            {{ SC_(50.0), SC_(0.125), SC_(0.25), SC_(0.196321043776719739372196642514913879497104766409504746018939) }},
            {{ SC_(1.125), SC_(1.0), SC_(0.25), SC_(1.77299767784815770192352979665283069318388205110727241629752) }},
            {{ SC_(1.125), SC_(10.0), SC_(0.25), SC_(0.662467818678976949597336360256848770217429434745967677192487) }},
            {{ SC_(1.125), SC_(3.0), SC_(0.25), SC_(-0.142697285116693775525461312178015106079842313950476205580178) }},
            {{ SC_(1.00390625) /*T(257)/256*/, SC_(1.5), SC_(0.125), SC_(22.2699300473528164111357290313578126108398839810535700884237) }},
            {{ SC_(1.00390625) /*T(257)/256*/, SC_(21.5), SC_(0.125), SC_(-0.535406081652313940727588125663856894154526187713506526799429) }},
            // Bug cases from Rocco Romeo:
            { { SC_(-1E-170), SC_(0.785398163397448309615660845819875721049292349843776455243) /*boost::math::constants::pi<T>() / 4*/, SC_(1E-164), SC_(0.785398163397448309615660845819875721049292349843776455243736) } },
            { { SC_(-1E-170), SC_(0.785398163397448309615660845819875721049292349843776455243) /*boost::math::constants::pi<T>() / 4*/, SC_(-1E-164), SC_(0.785398163397448309615660845819875721049292349843776455243736) } },
            { { SC_(-2.220446049250313080847263336181640625e-16) /*-ldexp(T(1.0), -52)*/, SC_(0.785398163397448309615660845819875721049292349843776455243) /*boost::math::constants::pi<T>() / 4*/, SC_(0.9375), SC_(0.866032844934895872810905364370384153285798081574191920571016) } },
            { { SC_(-2.220446049250313080847263336181640625e-16) /*-ldexp(T(1.0), -52)*/, SC_(0.785398163397448309615660845819875721049292349843776455243) /*boost::math::constants::pi<T>() / 4*/, SC_(-0.9375), SC_(0.866032844934895872810905364370384153285798081574191920571016) } },
            { { std::numeric_limits<T>::max_exponent > 600 ? SC_(-6.63922491020958873361985258190585784161991397158789903999305172750504449826065303423953127831536607086111661e181) /*-ldexp(T(1), 604)*/ : SC_(0.0),
                SC_(-2.2883557340936751629907904626893087059630763187242253081319190179713605588564915508590211920126569156532129e-246) /*-ldexp(T(1), -816)*/,
                SC_(2.9833362924800826973163861261851735349505886138401627577086073332825113540280473219744699561601711213446046e-154) /*ldexp(T(1), -510)*/, std::numeric_limits<T>::max_exponent > 600 ? SC_(-2.28835573409367516299079046268930870596307631872422530813192e-246) : SC_(-2.28835573409367516299079046268930870596307631872422530813192e-246) } },
            { { std::numeric_limits<T>::max_exponent > 600 ? SC_(-6.63922491020958873361985258190585784161991397158789903999305172750504449826065303423953127831536607086111661e181) /*-ldexp(T(1), 604)*/ : SC_(0.0),
                SC_(-2.2883557340936751629907904626893087059630763187242253081319190179713605588564915508590211920126569156532129e-246) /*-ldexp(T(1), -816)*/,
                SC_(-2.9833362924800826973163861261851735349505886138401627577086073332825113540280473219744699561601711213446046e-154) /*-ldexp(T(1), -510)*/,
                std::numeric_limits<T>::max_exponent > 600 ? SC_(-2.28835573409367516299079046268930870596307631872422530813192e-246) : SC_(-2.28835573409367516299079046268930870596307631872422530813192e-246) } },
            { { SC_(-5.7456966998645880645292998187840197954917072165248345657117366246194931336885953843013801519761445019343718e-188) /*-ldexp(T(1), -622)*/, SC_(-1.4996968138956309548176444376280653535399616962391082979373344476177108558521903027709681283974148362424896e-241) /*-ldexp(T(1), -800)*/, SC_(1.83670992315982423120115083940975887159166493245638675235742454106002696789801120758056640625e-40) /*ldexp(T(1), -132)*/, SC_(-1.49969681389563095481764443762806535353996169623910829793733e-241) } },
            { { SC_(-5.7456966998645880645292998187840197954917072165248345657117366246194931336885953843013801519761445019343718e-188) /*-ldexp(T(1), -622)*/, SC_(-1.4996968138956309548176444376280653535399616962391082979373344476177108558521903027709681283974148362424896e-241) /*-ldexp(T(1), -800)*/, SC_(-1.83670992315982423120115083940975887159166493245638675235742454106002696789801120758056640625e-40) /*-ldexp(T(1), -132)*/, SC_(-1.49969681389563095481764443762806535353996169623910829793733e-241) } },
            { { SC_(-6.6243372842224761350235742755846980122105019136217468782770353048106301097687150169478930287429355091446506e-170) /*-ldexp(T(1), -562)*/, SC_(7.1746481373430634031294954664443705921549411424077607513961896135157303433516062796115875244140625e-43) /*ldexp(T(1), -140)*/, SC_(8.6361685550944446253863518628003995711160003644362813850237034701685918031624270579715075034722882265605472e-78) /*ldexp(T(1), -256)*/, SC_(7.174648137343063403129495466444370592154941142407760751e-43) } },
            { { SC_(-6.6243372842224761350235742755846980122105019136217468782770353048106301097687150169478930287429355091446506e-170) /*-ldexp(T(1), -562)*/, SC_(-7.1746481373430634031294954664443705921549411424077607513961896135157303433516062796115875244140625e-43) /*-ldexp(T(1), -140)*/, SC_(8.6361685550944446253863518628003995711160003644362813850237034701685918031624270579715075034722882265605472e-78) /*ldexp(T(1), -256)*/, SC_(-7.17464813734306340312949546644437059215494114240776075e-43) } },
            { { SC_(7.7868710555449746371177365743005823355489124613059339568813918517897390050405261457702973319275391516778931e-208) /*ldexp(T(1), -688)*/, SC_(-7.0747492803333690371164994460060873286582274985462017106114178827621104051506602458902589468444985151984003e-74) /*-ldexp(T(1), -243)*/, SC_(7.9654595556622613851444019888385590279555227759630939303694292669308145075652908047154493736009134297049172e-59) /*ldexp(T(1), -193)*/, SC_(-7.07474928033336903711649944600608732865822749854620171e-74) } },
            { { SC_(-7.7868710555449746371177365743005823355489124613059339568813918517897390050405261457702973319275391516778931e-208) /*-ldexp(T(1), -688)*/, SC_(-7.0747492803333690371164994460060873286582274985462017106114178827621104051506602458902589468444985151984003e-74) /*-ldexp(T(1), -243)*/, SC_(7.9654595556622613851444019888385590279555227759630939303694292669308145075652908047154493736009134297049172e-59) /*ldexp(T(1), -193)*/, SC_(-7.07474928033336903711649944600608732865822749854620171e-74) } },
            // Special cases where k = 0:
            { { SC_(0.5), SC_(1.0), SC_(0.0), SC_(1.17881507892743738986863357869566288974084658835353613038547) } },
            { { SC_(-0.5), SC_(1.0), SC_(0.0), SC_(0.888286691263535380266337576823783210424994266596287990733270) } },
            { { SC_(0.5), SC_(-1.0), SC_(0.0), SC_(-1.17881507892743738986863357869566288974084658835353613038547) } },
            { { SC_(-0.5), SC_(-1.0), SC_(0.0), SC_(-0.888286691263535380266337576823783210424994266596287990733270) } },
            // k == 0 and phi > pi/2:
            { { SC_(0.5), SC_(2.0), SC_(0.0), SC_(3.03379730757207227653600089552126882582809860566558143254794) } },
            { { SC_(-0.5), SC_(2.0), SC_(0.0), SC_(1.57453655812023739911111328195028658229986230310938753315640) } },
            { { SC_(0.5), SC_(-2.0), SC_(0.0), SC_(-3.03379730757207227653600089552126882582809860566558143254794) } },
            { { SC_(-0.5), SC_(-2.0), SC_(0.0), SC_(-1.57453655812023739911111328195028658229986230310938753315640) } },
            // Special cases where k = 1:
            { { SC_(0.5), SC_(1.0), SC_(1.0), SC_(1.4830998734200773326887632776553375078936815318419194718912351) } },
            { { SC_(-0.5), SC_(1.0), SC_(1.0), SC_(1.07048347329000030842347009377117215811122412769516781788253) } },
            { { SC_(0.5), SC_(-1.0), SC_(1.0), SC_(-1.4830998734200773326887632776553375078936815318419194718912) } },
            { { SC_(-0.5), SC_(-1.0), SC_(1.0), SC_(-1.07048347329000030842347009377117215811122412769516781788253) } },
            // special cases where v = 1:
            { { SC_(1.0), SC_(0.5), SC_(0.5), SC_(0.55225234291197632914658859230278152249148960801635386133501) } },
            { { SC_(1.0), SC_(-0.5), SC_(0.5), SC_(-0.55225234291197632914658859230278152249148960801635386133501) } },
            { { SC_(1.0), SC_(2.0), SC_(0.5), SC_(-2.87534521505997989921579168327307068134740792740155171368532) } },
            { { SC_(1.0), SC_(-2.0), SC_(0.5), SC_(2.87534521505997989921579168327307068134740792740155171368532) } },
            { { SC_(1.0), SC_(2.0), SC_(6.2230152778611417071440640537801242405902521687211671331011166147896988340353834411839448231257136169569665e-61) /*ldexp(T(1), -200)*/, SC_(-2.18503986326151899164330610231368254343201774622766316456295) } },
            { { SC_(1.0), SC_(-2.0), SC_(6.2230152778611417071440640537801242405902521687211671331011166147896988340353834411839448231257136169569665e-61) /*ldexp(T(1), -200)*/, SC_(2.18503986326151899164330610231368254343201774622766316456295) } },
            { { SC_(1.0), SC_(7.00649232162408535461864791644958065640130970938257885878534141944895541342930300743319094181060791015625e-46) /*ldexp(T(1.0), -150)*/, SC_(6.2230152778611417071440640537801242405902521687211671331011166147896988340353834411839448231257136169569665e-61) /*ldexp(T(1), -200)*/, SC_(7.006492321624085354618647916449580656401309709382578858e-46) } },
            { { SC_(1.0), SC_(-7.00649232162408535461864791644958065640130970938257885878534141944895541342930300743319094181060791015625e-46) /*ldexp(T(1.0), -150)*/, SC_(-6.2230152778611417071440640537801242405902521687211671331011166147896988340353834411839448231257136169569665e-61) /*-ldexp(T(1), -200)*/, SC_(-7.006492321624085354618647916449580656401309709382578858e-46) } },
            // Previously unsupported region with v > 1 and |phi| > PI/2, this is the only region
            // with high-ish error rates caused by argument reduction by Pi:
            { { SC_(20.0), SC_(3.1425685882568359375) /*ldexp(T(1647611), -19)*/, SC_(0.5), SC_(0.000975940902692994840122139131147517258405256880370413541280) } },
            { { SC_(20.0), SC_(-3.1425685882568359375) /*-ldexp(T(1647611), -19)*/, SC_(0.5), SC_(-0.000975940902692994840122139131147517258405256880370413541280) } },
            { { SC_(1.015625) /*T(1.0) + ldexp(T(1), -6)*/, SC_(1.6957950592041015625) /*ldexp(T(889085), -19)*/, SC_(0.5), SC_(-27.1647225624906589308619292363045712770651414487085887109197) } },
            { { SC_(1.015625) /*T(1.0) + ldexp(T(1), -6)*/, SC_(-1.6957950592041015625) /*-ldexp(T(889085), -19)*/, SC_(0.5), SC_(27.1647225624906589308619292363045712770651414487085887109197) } },
            // Phi = 0:
            { { SC_(1.0), SC_(0.0), SC_(0.5), SC_(0.0) } },
            { { SC_(-1.0), SC_(0.0), SC_(0.5), SC_(0.0) } },
            { { SC_(100.0), SC_(0.0), SC_(0.5), SC_(0.0) } },
            { { SC_(-100.0), SC_(0.0), SC_(0.5), SC_(0.0) } },
        } };

#include "math/test/ellint_pi3_data.ipp"
#include "math/test/ellint_pi3_large_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[2], datum[0], datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[2], datum[0], datum[1]));
                if (!(actual == datum[3])) // +/-inf is equal to, but not "close" to, +/-inf
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[3], tolerance);
            };
        };

        ::for_each(data1, tester(600 * eps<T>));
        ::for_each(ellint_pi3_data, tester(9 * eps<T>));
        ::for_each(ellint_pi3_large_data, tester(3 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_ellint_3_boundaries, T, fptypes) {
        auto const tolerance = eps<T>;

        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1), static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        BOOST_CHECK(std::isnan(test_fn<T>(T(1.0001), T(-1), T(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(0.5), T(20), T(1.5))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(1), T(0.5), T(2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(1), T(-0.5), T(2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(1), T(-0.5), T(-2))));
        BOOST_CHECK(verify_domain_error());

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1), static_cast<T>(0), static_cast<T>(0)),
            static_cast<T>(0L), tolerance);
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(0), static_cast<T>(1), static_cast<T>(0)),
            static_cast<T>(0L), tolerance);
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace ellint_3

namespace expint {
    template<class T>
    constexpr auto control_fn = [](T x) {
        return boost::math::expint(x);
    };

    template<class T>
    constexpr auto test_fn = [](auto x) {
        static_assert(always_false<decltype(x)>);
    };
    template<>
    constexpr auto test_fn<float> = std::expintf;
    template<>
    constexpr auto test_fn<double> = std::expint;
    template<>
    constexpr auto test_fn<long double> = std::expintl;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_expint, T, fptypes) {
        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0]));
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[1], tolerance);
            };
        };

#include "math/test/expinti_data.ipp"
#include "math/test/expinti_data_double.ipp"

        ::for_each(expinti_data, tester(2 * eps<T>));

        if(boost::math::tools::log_max_value<T>() > 100) {
            ::for_each(expinti_data_double, tester(2 * eps<T>));
        }

        auto const tolerance = eps<T>;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1)/1024), static_cast<T>(-6.35327933972759151358547423727042905862963067106751711596065L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(0.125)), static_cast<T>(-1.37320852494298333781545045921206470808223543321810480716122L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(0.5)), static_cast<T>(0.454219904863173579920523812662802365281405554352642045162818L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(1)), static_cast<T>(1.89511781635593675546652093433163426901706058173270759164623L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(50.5)), static_cast<T>(1.72763195602911805201155668940185673806099654090456049881069e20L), tolerance);

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(-1)/1024), static_cast<T>(-6.35523246483107180261445551935803221293763008553775821607264L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(-0.125)), static_cast<T>(-1.62342564058416879145630692462440887363310605737209536579267L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(-0.5)), static_cast<T>(-0.559773594776160811746795939315085235226846890316353515248293L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(-1)), static_cast<T>(-0.219383934395520273677163775460121649031047293406908207577979L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(static_cast<T>(-50.5)), static_cast<T>(-2.27237132932219350440719707268817831250090574830769670186618e-24L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_expint_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        BOOST_CHECK_EQUAL(test_fn<T>(static_cast<T>(0)), -inf<T>);
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace expint

namespace hermite {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>);
    };
    template<>
    constexpr auto test_fn<float> = std::hermitef;
    template<>
    constexpr auto test_fn<double> = std::hermite;
    template<>
    constexpr auto test_fn<long double> = std::hermitel;

    template<class T>
    constexpr auto control_fn = [](unsigned n, T x) {
        return boost::math::hermite(n, x);
    };

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_hermite, T, fptypes) {
#include "math/test/hermite.ipp"

        for(auto const& datum : hermite) {
            unsigned const n = std::lround(datum[0]);
            auto const actual = test_fn<T>(n, datum[1]);
            BOOST_CHECK_EQUAL(actual, control_fn<T>(n, datum[1]));
            if (!(actual == datum[2])) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 5 * eps<T>);
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_hermite_spots, T, fptypes) {
        auto const tolerance = 2 * eps<T>;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(0, static_cast<T>(1)), static_cast<T>(1.L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(1)), static_cast<T>(2.L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(2)), static_cast<T>(4.L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(10)), static_cast<T>(20), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(100)), static_cast<T>(200), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(1e6)), static_cast<T>(2e6), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, static_cast<T>(30)), static_cast<T>(5.896624628001300E+17L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, static_cast<T>(1000)), static_cast<T>(1.023976960161280E+33L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, static_cast<T>(10)), static_cast<T>(8.093278209760000E+12L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(10, static_cast<T>(-10)), static_cast<T>(8.093278209760000E+12L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(3, static_cast<T>(-10)), static_cast<T>(-7.880000000000000E+3L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(3, static_cast<T>(-1000)), static_cast<T>(-7.999988000000000E+9L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(3, static_cast<T>(-1000000)), static_cast<T>(-7.999999999988000E+18L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_hermite_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(1u, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace hermite

namespace laguerre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>);
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
            if (!(actual == datum[2])) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 3200 * eps<T>);
        }
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_laguerre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(0.5L)),
            static_cast<T>(0.5L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(0.5L)),
            static_cast<T>(-0.3307291666666666666666666666666666666667L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(0.5L)),
            static_cast<T>(-0.5183392237103174603174603174603174603175L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(20, static_cast<T>(0.5L)),
            static_cast<T>(0.3120174870800154148915399248893113634676L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, static_cast<T>(0.5L)),
            static_cast<T>(-0.3181388060269979064951118308575628226834L), tolerance);

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(-0.5L)),
            static_cast<T>(1.5L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(-0.5L)),
            static_cast<T>(3.835937500000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(-0.5L)),
            static_cast<T>(7.950934709821428571428571428571428571429L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(20, static_cast<T>(-0.5L)),
            static_cast<T>(76.12915699869631476833699787070874048223L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, static_cast<T>(-0.5L)),
            static_cast<T>(2307.428631277506570629232863491518399720L), tolerance);

        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(4.5L)),
            static_cast<T>(-3.500000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(4.5L)),
            static_cast<T>(0.08593750000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(4.5L)),
            static_cast<T>(-1.036928013392857142857142857142857142857L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(20, static_cast<T>(4.5L)),
            static_cast<T>(1.437239150257817378525582974722170737587L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(50, static_cast<T>(4.5L)),
            static_cast<T>(-0.7795068145562651416494321484050019245248L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_laguerre_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(1u, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace laguerre

namespace legendre {
    template<class>
    constexpr auto test_fn = [](unsigned, auto x) {
        static_assert(always_false<decltype(x)>);
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
            if (!(actual == datum[2])) // +/-inf is equal to, but not "close" to, +/-inf
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], 310 * eps<T>);
        };

        ::for_each(legendre_p, test_one_value);
        ::for_each(legendre_p_large, test_one_value);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_legendre_spots, T, fptypes) {
        auto const tolerance = eps<T> * 100;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(1, static_cast<T>(0.5L)),
            static_cast<T>(0.5L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-1, static_cast<T>(0.5L)),
            static_cast<T>(1L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(4, static_cast<T>(0.5L)),
            static_cast<T>(-0.2890625000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-4, static_cast<T>(0.5L)),
            static_cast<T>(-0.4375000000000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(7, static_cast<T>(0.5L)),
            static_cast<T>(0.2231445312500000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-7, static_cast<T>(0.5L)),
            static_cast<T>(0.3232421875000000000000000000000000000000L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(40, static_cast<T>(0.5L)),
            static_cast<T>(-0.09542943523261546936538467572384923220258L), tolerance);
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(-40, static_cast<T>(0.5L)),
            static_cast<T>(-0.1316993126940266257030910566308990611306L), tolerance);
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_legendre_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(1u, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |x| <= 1
        BOOST_CHECK(std::isnan(test_fn<T>(0u, static_cast<T>(2))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(0u, static_cast<T>(-2))));
        BOOST_CHECK(verify_domain_error());

        auto const tolerance = eps<T>;
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(0u, static_cast<T>(1)), static_cast<T>(1), tolerance);
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_CLOSE_FRACTION(test_fn<T>(0u, static_cast<T>(-1)), static_cast<T>(1), tolerance);
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace legendre

namespace riemann_zeta {
    template<class T>
    constexpr auto control_fn = [](T x) {
        return boost::math::zeta(x);
    };

    template<class T>
    constexpr auto test_fn = [](auto x) {
        static_assert(always_false<decltype(x)>);
    };
    template<>
    constexpr auto test_fn<float> = std::riemann_zetaf;
    template<>
    constexpr auto test_fn<double> = std::riemann_zeta;
    template<>
    constexpr auto test_fn<long double> = std::riemann_zetal;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_riemann_zeta, T, fptypes) {
        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                auto const actual = test_fn<T>(datum[0]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(datum[0]));
                BOOST_CHECK_CLOSE_FRACTION(actual, datum[1], tolerance);
            };
        };

#include "math/test/zeta_data.ipp"
#include "math/test/zeta_neg_data.ipp"
#include "math/test/zeta_1_up_data.ipp"
#include "math/test/zeta_1_below_data.ipp"

        boost::array<boost::array<typename table_type<T>::type, 2>, 90> integer_data = { {
            {{ SC_(2.0), SC_(1.6449340668482264364724151666460252) }}, {{ SC_(3.0), SC_(1.2020569031595942853997381615114500) }}, {{ SC_(4.0), SC_(1.0823232337111381915160036965411679) }}, {{ SC_(5.0), SC_(1.0369277551433699263313654864570342) }}, {{ SC_(6.0), SC_(1.0173430619844491397145179297909205) }}, {{ SC_(7.0), SC_(1.0083492773819228268397975498497968) }}, {{ SC_(8.0), SC_(1.0040773561979443393786852385086525) }}, {{ SC_(9.0), SC_(1.0020083928260822144178527692324121) }}, {{ SC_(10.0), SC_(1.0009945751278180853371459589003190) }}, {{ SC_(11.0), SC_(1.0004941886041194645587022825264699) }}, {{ SC_(12.0), SC_(1.0002460865533080482986379980477397) }}, {{ SC_(13.0), SC_(1.0001227133475784891467518365263574) }}, {{ SC_(14.0), SC_(1.0000612481350587048292585451051353) }}, {{ SC_(15.0), SC_(1.0000305882363070204935517285106451) }}, {{ SC_(16.0), SC_(1.0000152822594086518717325714876367) }}, {{ SC_(17.0), SC_(1.0000076371976378997622736002935630) }}, {{ SC_(18.0), SC_(1.0000038172932649998398564616446219) }}, {{ SC_(19.0), SC_(1.0000019082127165539389256569577951) }}, {{ SC_(20.0), SC_(1.0000009539620338727961131520386834) }}, {{ SC_(21.0), SC_(1.0000004769329867878064631167196044) }}, {{ SC_(22.0), SC_(1.0000002384505027277329900036481868) }}, {{ SC_(23.0), SC_(1.0000001192199259653110730677887189) }}, {{ SC_(24.0), SC_(1.0000000596081890512594796124402079) }}, {{ SC_(25.0), SC_(1.0000000298035035146522801860637051) }}, {{ SC_(26.0), SC_(1.0000000149015548283650412346585066) }}, {{ SC_(27.0), SC_(1.0000000074507117898354294919810042) }}, {{ SC_(28.0), SC_(1.0000000037253340247884570548192040) }}, {{ SC_(29.0), SC_(1.0000000018626597235130490064039099) }}, {{ SC_(30.0), SC_(1.0000000009313274324196681828717647) }}, {{ SC_(31.0), SC_(1.0000000004656629065033784072989233) }}, {{ SC_(32.0), SC_(1.0000000002328311833676505492001456) }}, {{ SC_(33.0), SC_(1.0000000001164155017270051977592974) }}, {{ SC_(34.0), SC_(1.0000000000582077208790270088924369) }}, {{ SC_(35.0), SC_(1.0000000000291038504449709968692943) }}, {{ SC_(36.0), SC_(1.0000000000145519218910419842359296) }}, {{ SC_(37.0), SC_(1.0000000000072759598350574810145209) }}, {{ SC_(38.0), SC_(1.0000000000036379795473786511902372) }}, {{ SC_(39.0), SC_(1.0000000000018189896503070659475848) }}, {{ SC_(40.0), SC_(1.0000000000009094947840263889282533) }}, {{ SC_(41.0), SC_(1.0000000000004547473783042154026799) }}, {{ SC_(42.0), SC_(1.0000000000002273736845824652515227) }}, {{ SC_(43.0), SC_(1.0000000000001136868407680227849349) }}, {{ SC_(44.0), SC_(1.0000000000000568434198762758560928) }}, {{ SC_(45.0), SC_(1.0000000000000284217097688930185546) }}, {{ SC_(46.0), SC_(1.0000000000000142108548280316067698) }}, {{ SC_(47.0), SC_(1.0000000000000071054273952108527129) }}, {{ SC_(48.0), SC_(1.0000000000000035527136913371136733) }}, {{ SC_(49.0), SC_(1.0000000000000017763568435791203275) }}, {{ SC_(50.0), SC_(1.0000000000000008881784210930815903) }}, {{ SC_(51.0), SC_(1.0000000000000004440892103143813364) }}, {{ SC_(52.0), SC_(1.0000000000000002220446050798041984) }}, {{ SC_(53.0), SC_(1.0000000000000001110223025141066134) }}, {{ SC_(54.0), SC_(1.0000000000000000555111512484548124) }}, {{ SC_(55.0), SC_(1.0000000000000000277555756213612417) }}, {{ SC_(56.0), SC_(1.0000000000000000138777878097252328) }}, {{ SC_(57.0), SC_(1.0000000000000000069388939045441537) }}, {{ SC_(58.0), SC_(1.0000000000000000034694469521659226) }}, {{ SC_(59.0), SC_(1.0000000000000000017347234760475766) }}, {{ SC_(60.0), SC_(1.0000000000000000008673617380119934) }},
            {{ SC_(-61.0), SC_(-3.3066089876577576725680214670439210e34) }}, {{ SC_(-59.0), SC_(3.5666582095375556109684574608651829e32) }}, {{ SC_(-57.0), SC_(-4.1147288792557978697665486067619336e30) }}, {{ SC_(-55.0), SC_(5.0890659468662289689766332915911925e28) }}, {{ SC_(-53.0), SC_(-6.7645882379292820990945242301798478e26) }}, {{ SC_(-51.0), SC_(9.6899578874635940656497942894654088e24) }}, {{ SC_(-49.0), SC_(-1.5001733492153928733711440151515152e23) }}, {{ SC_(-47.0), SC_(2.5180471921451095697089023320225526e21) }}, {{ SC_(-45.0), SC_(-4.5979888343656503490437943262411348e19) }}, {{ SC_(-43.0), SC_(9.1677436031953307756992753623188406e17) }}, {{ SC_(-41.0), SC_(-2.0040310656516252738108421663238939e16) }}, {{ SC_(-39.0), SC_(4.8241448354850170371581670362158167e14) }}, {{ SC_(-37.0), SC_(-1.2850850499305083333333333333333333e13) }}, {{ SC_(-35.0), SC_(3.8087931125245368811553022079337869e11) }}, {{ SC_(-33.0), SC_(-1.2635724795916666666666666666666667e10) }}, {{ SC_(-31.0), SC_(4.7238486772162990196078431372549020e8) }}, {{ SC_(-29.0), SC_(-2.0052695796688078946143462272494531e7) }}, {{ SC_(-27.0), SC_(974936.82385057471264367816091954023) }}, {{ SC_(-25.0), SC_(-54827.583333333333333333333333333333) }}, {{ SC_(-23.0), SC_(3607.5105463980463980463980463980464) }}, {{ SC_(-21.0), SC_(-281.46014492753623188405797101449275) }}, {{ SC_(-19.0), SC_(26.456212121212121212121212121212121) }}, {{ SC_(-17.0), SC_(-3.0539543302701197438039543302701197) }}, {{ SC_(-15.0), SC_(0.44325980392156862745098039215686275) }}, {{ SC_(-13.0), SC_(-0.083333333333333333333333333333333333) }}, {{ SC_(-11.0), SC_(0.021092796092796092796092796092796093) }}, {{ SC_(-9.0), SC_(-0.0075757575757575757575757575757575758) }}, {{ SC_(-7.0), SC_(0.0041666666666666666666666666666666667) }}, {{ SC_(-5.0), SC_(-0.0039682539682539682539682539682539683) }}, {{ SC_(-3.0), SC_(0.0083333333333333333333333333333333333) }}, {{ SC_(-1.0), SC_(-0.083333333333333333333333333333333333) }}
        } };

        ::for_each(zeta_data, tester(eps<T>));
        ::for_each(zeta_neg_data, tester(10 * eps<T>));
        ::for_each(zeta_1_up_data, tester(eps<T>));
        ::for_each(zeta_1_below_data, tester(eps<T>));
        ::for_each(integer_data, tester(10 * eps<T>));

        //
        // Basic sanity checks, tolerance is either 5 or 10 epsilon
        // expressed as a percentage:
        //
        BOOST_MATH_STD_USING
        T tolerance = boost::math::tools::epsilon<T>() * 100 *
            (boost::is_floating_point<T>::value ? 5 : 10);
        // An extra fudge factor for real_concept which has a less accurate tgamma:
        T tolerance_tgamma_extra = std::numeric_limits<T>::is_specialized ? 1 : 10;

        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(0.125)), static_cast<T>(-0.63277562349869525529352526763564627152686379131122L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(1023) / static_cast<T>(1024)), static_cast<T>(-1023.4228554489429786541032870895167448906103303056L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(1025) / static_cast<T>(1024)), static_cast<T>(1024.5772867695045940578681624248887776501597556226L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(0.5)), static_cast<T>(-1.46035450880958681288949915251529801246722933101258149054289L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(1.125)), static_cast<T>(8.5862412945105752999607544082693023591996301183069L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(2)), static_cast<T>(1.6449340668482264364724151666460251892189499012068L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(3.5)), static_cast<T>(1.1267338673170566464278124918549842722219969574036L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(4)), static_cast<T>(1.08232323371113819151600369654116790277475095191872690768298L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(4 + static_cast<T>(1) / 1024), static_cast<T>(1.08225596856391369799036835439238249195298434901488518878804L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(4.5)), static_cast<T>(1.05470751076145426402296728896028011727249383295625173068468L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(6.5)), static_cast<T>(1.01200589988852479610078491680478352908773213619144808841031L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(7.5)), static_cast<T>(1.00582672753652280770224164440459408011782510096320822989663L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(8.125)), static_cast<T>(1.0037305205308161603183307711439385250181080293472L), 2 * tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(16.125)), static_cast<T>(1.0000140128224754088474783648500235958510030511915L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(0)), static_cast<T>(-0.5L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-0.125)), static_cast<T>(-0.39906966894504503550986928301421235400280637468895L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-1)), static_cast<T>(-0.083333333333333333333333333333333333333333333333333L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-2)), static_cast<T>(0L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-2.5)), static_cast<T>(0.0085169287778503305423585670283444869362759902200745L), tolerance * 3);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-3)), static_cast<T>(0.0083333333333333333333333333333333333333333333333333L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-4)), static_cast<T>(0), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-20)), static_cast<T>(0), tolerance * 100);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-21)), static_cast<T>(-281.46014492753623188405797101449275362318840579710L), tolerance * 100);
        BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-30.125)), static_cast<T>(2.2762941726834511267740045451463455513839970804578e7L), tolerance * 100);
        // Very small values:
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -20)), static_cast<T>(-0.500000876368989859479646132126454890645615288202492097957612L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -21)), static_cast<T>(-0.500000438184266833093492063649184012943132422189989164545507L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -22)), static_cast<T>(-0.500000219092076392425852854644256723571669269957526445270374L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -23)), static_cast<T>(-0.500000109546023940187789325464529558825433290921168958481804L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -24)), static_cast<T>(-0.500000054773008406088246161057525197302821575823476487961574L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -25)), static_cast<T>(-0.500000027386503312042790426817221131071450407798601059264341L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -26)), static_cast<T>(-0.500000013693251433271071983943082871935521396740331377486886L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -27)), static_cast<T>(-0.500000006846625660947956426350389518286874288247134329498289L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -28)), static_cast<T>(-0.500000003423312816552083476988056486473169377162409806781384L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -29)), static_cast<T>(-0.500000001711656404795568073849512135664960180586820144333542L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -30)), static_cast<T>(-0.500000000855828201527665623188910582717329375986726355164261L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -31)), static_cast<T>(-0.500000000427914100546303208463654361814800355929815322493143L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -32)), static_cast<T>(-0.500000000213957050218769203487022003676593508474107873788445L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -33)), static_cast<T>(-0.500000000106978525095789001562046589421133388262409441738089L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(ldexp(static_cast<T>(1), -34)), static_cast<T>(-0.500000000053489262544495600736249301842352101231724731340202L), tolerance);

        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -20)), static_cast<T>(-0.499999123632834911086872289657767335473025908373776645822722L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -21)), static_cast<T>(-0.499999561816189359548137231641582253243376087534976981434190L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -22)), static_cast<T>(-0.499999780908037655734554449793729262345041281451929584703788L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -23)), static_cast<T>(-0.499999890454004571852312499433422838864632848598847415933664L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -24)), static_cast<T>(-0.499999945226998721921779295091241395945379526155584220813497L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -25)), static_cast<T>(-0.499999972613498469959715937215237923104705216198368099221577L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -26)), static_cast<T>(-0.499999986306749012229554607064736104475024094525587925697276L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -27)), static_cast<T>(-0.499999993153374450427200221401546739119918746163907954406855L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -28)), static_cast<T>(-0.499999996576687211291705684949926422460038672790251466963619L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -29)), static_cast<T>(-0.499999998288343602165379216634983519354686193860717726606017L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -30)), static_cast<T>(-0.499999999144171800212571199432213326524228740247618955829902L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -31)), static_cast<T>(-0.499999999572085899888755997191626615213504580792674808876724L), tolerance * tolerance_tgamma_extra);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -32)), static_cast<T>(-0.499999999786042949889995597926798240562852438685508646794693L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -33)), static_cast<T>(-0.499999999893021474931402198791408471637626205588681812641711L), tolerance);
        BOOST_CHECK_CLOSE(::boost::math::zeta(-ldexp(static_cast<T>(1), -34)), static_cast<T>(-0.499999999946510737462302199352114463422268928922372277519378L), tolerance);
#ifdef BOOST_MSVC
#pragma warning(push)
#pragma warning(disable:4127 4756)
#endif
        //
        // Very large negative values need special handling in our code, test them here, due to a bug report by Rocco Romeo:
        //
        BOOST_CHECK_EQUAL(::boost::math::zeta(static_cast<T>(-200)), static_cast<T>(0));
        if(std::numeric_limits<T>::max_exponent >= 1024)
        {
            BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-171)), static_cast<T>(1.28194898634822427378088228065956967928127061276520385040051e172L), tolerance * 200);
            BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-171.5)), static_cast<T>(4.73930233055054501360661283732419615206017226423071857829425e172L), tolerance * 1000);
            BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-172.5)), static_cast<T>(-1.30113885243175165293156588942160456456090687128236657847674e174L), tolerance * 100);
            BOOST_CHECK_CLOSE(::boost::math::zeta(static_cast<T>(-173)), static_cast<T>(-9.66241211085609184243169684777934860657838245104636064505158e174L), tolerance * 100);
        }
        if(std::numeric_limits<T>::has_infinity)
        {
            BOOST_CHECK_EQUAL(::boost::math::zeta(static_cast<T>(-10007)), std::numeric_limits<T>::infinity());
            BOOST_CHECK_EQUAL(::boost::math::zeta(static_cast<T>(-10009)), -std::numeric_limits<T>::infinity());
        }
#ifdef BOOST_MSVC
#pragma warning(pop)
#endif
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_riemann_zeta_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());
    }
} // namespace riemann_zeta

namespace sph_bessel {
    template<class T>
    constexpr auto control_fn = [](unsigned n, T x) {
        return boost::math::sph_bessel(n, x);
    };

    template<class T>
    constexpr auto test_fn = [](auto n, T) {
        static_assert(always_false<decltype(n)>);
    };
    template<>
    constexpr auto test_fn<float> = std::sph_besself;
    template<>
    constexpr auto test_fn<double> = std::sph_bessel;
    template<>
    constexpr auto test_fn<long double> = std::sph_bessell;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_j, T, fptypes) {
#include "math/test/sph_bessel_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                unsigned const n = std::lround(datum[0]);
                auto const actual = test_fn<T>(n, datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(n, datum[1]));
                if (!(actual == datum[2])) {
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
                }
            };
        };

        ::for_each(sph_bessel_data, tester(250 * eps<T>));

        //
        // Some special cases:
        //
        BOOST_CHECK_EQUAL(test_fn<T>(0, T(0)), T(1));
        BOOST_CHECK_EQUAL(test_fn<T>(1, T(0)), T(0));
        BOOST_CHECK_EQUAL(test_fn<T>(100000, T(0)), T(0));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_cyl_bessel_j_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(1u, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        // domain is x >= 0
        BOOST_CHECK(!std::isnan(test_fn<T>(1u, T(0))));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(1u, -eps<T>)));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace sph_bessel

namespace sph_neumann {
    template<class T>
    constexpr auto control_fn = [](unsigned n, T x) {
        return boost::math::sph_neumann(n, x);
    };

    template<class T>
    constexpr auto test_fn = [](auto n, T) {
        static_assert(always_false<decltype(n)>);
    };
    template<>
    constexpr auto test_fn<float> = std::sph_neumannf;
    template<>
    constexpr auto test_fn<double> = std::sph_neumann;
    template<>
    constexpr auto test_fn<long double> = std::sph_neumannl;

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_sph_neumann, T, fptypes) {
#include "math/test/sph_neumann_data.ipp"

        auto const tester = [](T tolerance) {
            return [tolerance](auto const& datum) {
                unsigned const n = std::lround(datum[0]);
                auto const actual = test_fn<T>(n, datum[1]);
                BOOST_CHECK_EQUAL(actual, control_fn<T>(n, datum[1]));
                if (!(actual == datum[2])) {
                    BOOST_CHECK_CLOSE_FRACTION(actual, datum[2], tolerance);
                }
            };
        };

        ::for_each(sph_neumann_data, tester(290 * eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_sph_neumann_boundaries, T, fptypes) {
        errno = 0;
        BOOST_CHECK(std::isnan(test_fn<T>(1u, qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

        // domain is x >= 0

        //
        // Special cases that are errors:
        //
        BOOST_CHECK(!std::isnan(test_fn<T>(1u, T(0))));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(1u, -eps<T>)));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace sph_neumann

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
