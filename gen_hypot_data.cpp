#include <cassert>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <random>
#include <string>
#include <utility>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/tools/precision.hpp>

namespace {
    using big_float = boost::multiprecision::number<boost::multiprecision::backends::cpp_bin_float<1000,
        boost::multiprecision::backends::digit_base_2>>;

    big_float hypot3(double x, double y, double z) {
        x = std::abs(x);
        assert(!std::isinf(x) && !std::isnan(x));

        y = std::abs(y);
        assert(!std::isinf(y) && !std::isnan(y));

        z = std::abs(z);
        assert(!std::isinf(z) && !std::isnan(z));

        if (y > x) {
            std::swap(x, y);
        }

        if (z > x) {
            std::swap(x, z);
        }

        auto const frac_sq = [](double const num, double const denom) {
            big_float result = num;
            result /= denom;
            result *= result;
            return result;
        };

        return x * sqrt(1 + frac_sq(y, x) + frac_sq(z, x));
    }

    void generate(double const low, char const * const name) {
        constexpr auto n = 500;
        static auto engine = std::mt19937{};

        auto in_range = [low, high = low / boost::math::tools::epsilon<double>()](double const d) {
            return low <= d && d < high;
        };

        auto dist = std::uniform_int_distribution<unsigned long long>{};
        auto gen = [&] {
            double d;
            // Choose uniformly from the set of double representations in [low, high)
            do {
                auto const ull = dist(engine);
                std::memcpy(&d, &ull, sizeof(ull));
            } while(!in_range(std::fabs(d)));
            return d;
        };

        std::string basename = "hypot_" + std::string{name} + "_data";
        std::ofstream f{basename + ".ipp"};
        assert(f);
        f.precision(std::numeric_limits<double>::digits10);
        f.flags(std::ios_base::fmtflags(std::ios_base::scientific));
        f << "static const boost::array<boost::array<typename table_type<T>::type, 4>, " << n << "> " << basename << " = {{\n";

        for (auto i = 0; i < n; ++i) {
            double const d[3] = {gen(), gen(), gen()};
            f << "    {{ SC_(" << d[0] << "), SC_(" << d[1] << "), SC_(" << d[2] << "), SC_("
              << hypot3(d[0], d[1], d[2]) << ") }},\n";
        }

        f << "}};\n";
    }
}

int main() {
    generate(1.0, "low");
    generate(1.0e64, "high");
}
