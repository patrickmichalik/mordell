// Independently find integral points within a bounded coordinate range.
// Usage: box LO HI XMAX OUT.tsv, with -10^7 <= LO <= HI <= 10^7
// and 0 <= XMAX <= 10^8.

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

using Int128 = __int128_t;
using Int64 = long long;

Int64 root(Int128 n) {
    if (n < 0)
        throw std::runtime_error("Cannot take the square root of a negative integer.");

    Int64 r = (Int64)std::sqrt((long double)n);
    while (Int128(r) * r > n)
        --r;
    while (Int128(r + 1) * (r + 1) <= n)
        ++r;

    return r;
}

int main(int argc, char **argv) {
    try {
        if (argc != 5)
            throw std::runtime_error("Usage: box LO HI XMAX OUT.tsv.");
        Int64 lo = std::stoll(argv[1]), hi = std::stoll(argv[2]), xm = std::stoll(argv[3]);
        if (lo < -10000000 || lo > hi || hi > 10000000 || xm < 0 || xm > 100000000)
            throw std::runtime_error(
                "Bounds must satisfy -10^7 <= LO <= HI <= 10^7 and 0 <= XMAX <= 10^8.");

        std::ofstream out(argv[4]);
        if (!out)
            throw std::runtime_error("Could not open the output file.");

        Int64 xmin = 0;
        while (Int128(xmin - 1) * (xmin - 1) * (xmin - 1) + hi >= 0)
            --xmin;

        Int64 n = 0;
        for (Int64 x = xmin; x <= xm; ++x) {
            Int128 cube = Int128(x) * x * x, low = cube + lo, high = cube + hi;
            if (high < 0)
                continue;
            if (low < 0)
                low = 0;

            Int64 first = root(low);
            if (Int128(first) * first < low)
                ++first;
            Int64 last = root(high);

            for (Int64 y = first; y <= last; ++y) {
                Int64 k = (Int64)(Int128(y) * y - cube);
                if (k == 0)
                    continue;

                out << k << '\t' << x << '\t' << y << '\n';
                ++n;
            }
        }

        out.flush();
        if (!out)
            throw std::runtime_error("Could not write the output file.");

        std::cout << "{\"lower_inclusive\":" << lo << ",\"upper_inclusive\":" << hi
                  << ",\"x_min\":" << xmin << ",\"x_max\":" << xm
                  << ",\"points_y_nonnegative\":" << n << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
