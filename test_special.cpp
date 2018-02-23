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

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
