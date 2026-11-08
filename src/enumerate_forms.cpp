// Exhaustive interval enumeration using exact integer arithmetic.
// F(u, v) = a * u^3 + 3 * b * u^2 * v + 3 * c * u * v^2 + d * v^3.
// The discriminant is disc(F) = -108 * k.
// Use the conservative converted Belabas bounds from the independent solver.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

using Int128 = __int128_t;
using Int64 = long long;

Int64 floor_div(Int128 a, Int64 b) {
    if (b <= 0)
        throw std::runtime_error("The denominator must be positive.");

    Int128 q = a / b, r = a % b;
    return (Int64)(q - (r < 0));
}

Int64 ceil_div(Int128 a, Int64 b) { return -floor_div(-a, b); }

Int64 isqrt(Int128 n) {
    if (n < 0)
        throw std::runtime_error("Cannot take the square root of a negative integer.");

    Int64 r = (Int64)std::sqrt((long double)n);
    while (Int128(r + 1) * (r + 1) <= n)
        ++r;
    while (Int128(r) * r > n)
        --r;

    return r;
}

Int64 ceil_cbrt(Int128 n) {
    if (n < 0)
        throw std::runtime_error("Cannot take the cube root of a negative integer.");

    Int64 r = 0, h = 1;
    while (Int128(h) * h * h < n)
        h *= 2;

    while (r < h) {
        Int64 m = r + (h - r) / 2;
        if (Int128(m) * m * m >= n)
            h = m;
        else
            r = m + 1;
    }

    return r;
}

struct Counts {
    Int64 abc = 0, d = 0, forms = 0, reducible = 0;
};

void enumerate(Int64 lo, Int64 hi, std::ofstream &forms, Counts &counts) {
    if (lo > hi)
        return;

    bool positive_discriminant = hi < 0;
    Int64 T = 108 * std::max(std::abs(lo), std::abs(hi));
    Int64 A = 0, denom = positive_discriminant ? 729 : 27;
    while (Int128(A + 1) * (A + 1) * (A + 1) * (A + 1) * denom <= Int128(16) * T)
        ++A;

    for (Int64 a = 1; a <= A; ++a) {
        Int64 S = isqrt(positive_discriminant ? T : T / 3) + 1;
        Int64 b_max = (a + 1) / 2 + (isqrt(S) + 3) / 3 + 1;

        for (Int64 b = 0; b <= b_max; ++b) {
            Int64 c_min, c_max;
            if (positive_discriminant) {
                c_min = ceil_div(Int128(b) * b - isqrt(T / 81), a);
                c_max = b - a;
            } else {
                c_min = 1 - b;
                Int64 W = ceil_cbrt(ceil_div(T, 4 * a));
                c_max = a >= 2 * b ? floor_div(Int128(W) * a + 3 * Int128(b) * b, 3 * a)
                                   : floor_div(4 * Int128(W) + 12 * b - 3 * a, 12);
            }

            for (Int64 c = c_min; c <= c_max; ++c) {
                ++counts.abc;
                Int128 aa = Int128(a) * a, bb = 4 * Int128(b) * b * b - 6 * Int128(a) * b * c,
                       cc = 4 * Int128(a) * c * c * c - 3 * Int128(b) * b * c * c;
                Int128 delta = bb * bb - 4 * aa * (cc - 4 * Int128(hi));
                if (delta < 0)
                    continue;

                Int64 root = isqrt(delta);
                Int64 d_min = ceil_div(-Int128(root) - bb, 2 * a * a),
                      d_max = floor_div(Int128(root) - bb, 2 * a * a);

                Int64 P = b * b - a * c;
                if (positive_discriminant) {
                    if (P <= 0)
                        continue;
                    d_min = std::max(d_min, ceil_div(Int128(b) * c - P, a));
                    d_max = std::min(d_max, floor_div(Int128(b) * c + P, a));
                    if (b > 0)
                        d_max = std::min(d_max, floor_div(Int128(c) * c - P, b));
                } else {
                    Int64 low = -(a - 3 * b) * (a - 3 * b) - 3 * a * c,
                          up = (a + 3 * b) * (a + 3 * b) + 3 * a * c;
                    d_min = std::max(d_min, floor_div(Int128(low) + 9 * Int128(b) * c, a) + 1);
                    d_max = std::min(d_max, ceil_div(Int128(up) + 9 * Int128(b) * c, a) - 1);
                    if (b == 0)
                        d_min = std::max(d_min, 1LL);
                }

                for (Int64 d = d_min; d <= d_max; ++d) {
                    ++counts.d;
                    Int128 q = aa * d * d + bb * d + cc;
                    if (q % 4)
                        continue;
                    Int128 k = q / 4;
                    if (k < lo || k > hi || k == 0)
                        continue;
                    if (std::gcd(std::gcd(a, 3 * b), std::gcd(3 * c, d)) != 1)
                        continue;

                    if (positive_discriminant) {
                        Int64 Q = b * c - a * d, R = c * c - b * d;
                        if (!(std::abs(Q) <= P && P <= R))
                            continue;
                    } else {
                        if (!(Int128(d) * d - a * a > 3 * (Int128(b) * d - Int128(a) * c)))
                            continue;
                        Int128 middle = Int128(a) * d - 9 * Int128(b) * c;
                        if (!(-Int128(a - 3 * b) * (a - 3 * b) - 3 * Int128(a) * c < middle &&
                              middle < Int128(a + 3 * b) * (a + 3 * b) + 3 * Int128(a) * c))
                            continue;
                    }

                    forms << (Int64)k << '\t' << a << '\t' << b << '\t' << c << '\t' << d << '\n';
                    ++counts.forms;
                }
            }
        }
    }
}

int main(int argc, char **argv) {
    try {
        if (argc != 5)
            throw std::runtime_error(
                "Usage: enumerate LO HI FORMS.tsv REDUCIBLE.tsv (inclusive bounds).");
        Int64 lo = std::stoll(argv[1]), hi = std::stoll(argv[2]);
        if (lo > hi || lo < -10000000 || hi > 10000000)
            throw std::runtime_error("The range must satisfy -10^7 <= LO <= HI <= 10^7.");

        std::ofstream fs(argv[3]), rs(argv[4]);
        if (!fs || !rs)
            throw std::runtime_error("Could not open the output files.");

        Counts counts;
        enumerate(lo, std::min(hi, -1LL), fs, counts);
        enumerate(std::max(lo, 1LL), hi, fs, counts);

        Int64 c_max = isqrt(4 * Int128(std::max(std::abs(lo), std::abs(hi))));
        for (Int64 C = -c_max; C <= c_max; ++C) {
            if (!C)
                continue;

            Int128 num = 4 * Int128(C) * C * C - 4 * Int128(lo);
            if (num < 0)
                continue;

            Int64 b_max = isqrt(num / (3 * Int128(C) * C));
            for (Int64 B = 0; B <= b_max; ++B) {
                if ((B * C) % 2)
                    continue;

                Int128 numerator = Int128(C) * C * (4 * C - 3 * Int128(B) * B);
                if (numerator % 4)
                    throw std::runtime_error("The reducible-form numerator is not divisible by 4.");
                Int128 k = numerator / 4;
                if (k < lo || k > hi || k == 0)
                    continue;

                Int128 x = Int128(B) * B - C, y = Int128(B) * (3 * C - 2 * Int128(B) * B) / 2;
                if (y < 0)
                    y = -y;
                if (y * y - x * x * x != k)
                    throw std::runtime_error("The reducible-form point does not "
                                             "satisfy the curve equation.");

                rs << (Int64)k << '\t' << (Int64)x << '\t' << (Int64)y << '\t' << B << '\t' << C
                   << '\n';
                ++counts.reducible;
            }
        }

        fs.flush();
        rs.flush();
        if (!fs || !rs)
            throw std::runtime_error("Could not write the output files.");

        std::cout << "{\"lower_inclusive\":" << lo << ",\"upper_inclusive\":" << hi
                  << ",\"abc_candidates\":" << counts.abc << ",\"d_candidates\":" << counts.d
                  << ",\"candidate_forms\":" << counts.forms
                  << ",\"reducible_records\":" << counts.reducible << "}\n";
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
