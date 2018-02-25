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
