#define BOOST_TEST_MODULE WIP
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
        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(0))));
        BOOST_CHECK(verify_not_domain_error());

        // domain is |k| <= 1
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(2), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(-2), static_cast<T>(0))));
        BOOST_CHECK(verify_domain_error());

        BOOST_CHECK_EQUAL(test_fn<T>(static_cast<T>(1), static_cast<T>(0)), static_cast<T>(1));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK_EQUAL(test_fn<T>(static_cast<T>(-1), static_cast<T>(0)), static_cast<T>(1));
        BOOST_CHECK(verify_not_domain_error());

        BOOST_CHECK(std::isnan(test_fn<T>(T(1.0001), T(-1))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(0.5), T(1))));
        BOOST_CHECK(verify_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(T(0.5), T(2))));
        BOOST_CHECK(verify_domain_error());
    }
} // namespace comp_ellint_3

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

        ::for_each(data1, tester(eps<T>));
        // FIXME: ::for_each(ellint_pi3_data, tester(eps<T>));
        // FIXME: ::for_each(ellint_pi3_large_data, tester(eps<T>));
    }

    BOOST_AUTO_TEST_CASE_TEMPLATE(test_ellint_3_boundaries, T, fptypes) {
        auto const tolerance = eps<T>;

        errno = 0;
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

        BOOST_CHECK(std::isnan(test_fn<T>(qNaN<T>, static_cast<T>(1), static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), qNaN<T>, static_cast<T>(1))));
        BOOST_CHECK(verify_not_domain_error());
        BOOST_CHECK(std::isnan(test_fn<T>(static_cast<T>(1), static_cast<T>(1), qNaN<T>)));
        BOOST_CHECK(verify_not_domain_error());

#if 0 // FIXME
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
#endif
    }
} // namespace ellint_3

int main(int argc, char *argv[]) {
    auto const result = boost::unit_test::unit_test_main(init_unit_test, argc, argv);
    return result == boost::exit_success ? PM_TEST_PASS : PM_TEST_FAIL;
}
