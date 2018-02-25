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

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
