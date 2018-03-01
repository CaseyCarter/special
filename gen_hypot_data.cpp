#define _CRT_SECURE_NO_WARNINGS

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

    big_float hypot3(big_float x, big_float y, big_float z)
        {
        using std::swap;

        x = abs(x);
        assert(!isinf(x) && !isnan(x));

        y = abs(y);
        assert(!isinf(y) && !isnan(y));

        z = abs(z);
        assert(!isinf(z) && !isnan(z));

        if (y > x)
            {
            swap(x, y);
            }

        if (z > x)
            {
            swap(x, z);
            }

        const auto frac_sq = [](const big_float num, const big_float denom)
            {
            const big_float result = num / denom;
            return result * result;
            };

        return x * sqrt(1 + frac_sq(y, x) + frac_sq(z, x));
        }

    void generate(double const low, std::string name) {
        constexpr auto n = 500;
        static auto engine = std::mt19937{};

        double const high = low / boost::math::tools::epsilon<double>();

        auto dist = std::uniform_int_distribution<unsigned long long>{};
        auto gen = [&] {
            auto in_range = [=](double const d) {
                return low <= d && d < high;
            };
            double d;
            do {
                auto const ull = dist(engine);
                std::memcpy(&d, &ull, sizeof(ull));
            } while(!in_range(std::fabs(d)));
            return d;
        };

        name = "hypot_" + name + "_data";
        std::ofstream f{name + ".ipp"};
        assert(f);
        f.precision(std::numeric_limits<double>::digits10);
        f.flags(std::ios_base::fmtflags(std::ios_base::scientific));
        f << "static const boost::array<boost::array<typename table_type<T>::type, 4>, " << n << "> " << name << " = {{\n";

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
