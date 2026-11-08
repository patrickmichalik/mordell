#!/usr/bin/env python3
"""Exact counterexamples to the printed Bennett–Ghadermarzi reduction formulas.

Coefficient convention: F(u, v) = a * u^3 + 3 * b * u^2 * v + 3 * c * u * v^2 + d * v^3.
Print exact counterexamples as JSON using the Python standard library.
"""

import json
from fractions import Fraction

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


def discriminant(f):
    a, b, c, d = f
    return -27 * (
        a * a * d * d - 6 * a * b * c * d - 3 * b * b * c * c + 4 * a * c**3 + 4 * b**3 * d
    )


def is_reduced_positive(f, printed):
    a, b, c, d = f
    P = b * b - a * c
    Q = b * c - a * d
    R = c * c - b * d
    return (
        discriminant(f) > 0
        and abs(Q) <= P <= R
        and a > 0
        and b >= 0
        and (b != 0 or d < 0)
        and (Q != 0 or d < 0)
        and (P != Q or b < abs(a - b))
        and (P != R or (a <= abs(d) and (b < abs(c) if printed or a == abs(d) else True)))
    )


def is_reduced_negative(f, printed):
    a, b, c, d = f
    middle = 3 * (a * d - b * c) if printed else a * d - 9 * b * c
    return (
        discriminant(f) < 0
        and d * d - a * a > 3 * (b * d - a * c)
        and -((a - 3 * b) ** 2) - 3 * a * c < middle < (a + 3 * b) ** 2 + 3 * a * c
        and a > 0
        and b >= 0
        and (b != 0 or d > 0)
    )


def shift_variable(f, n):
    a, b, c, d = f
    return (
        a,
        b + a * n,
        c + 2 * b * n + a * n * n,
        d + 3 * c * n + 3 * b * n * n + a * n**3,
    )


def is_monic_irreducible(f):
    a, b, c, d = f
    assert a == 1 and d
    for z in range(-abs(d), abs(d) + 1):
        if z and d % z == 0 and a * z**3 + 3 * b * z * z + 3 * c * z + d == 0:
            return False
    return True


def main():
    d_example = (1, 3, 1, -3)
    assert (
        discriminant(d_example) == 7668
        and is_monic_irreducible(d_example)
        and is_reduced_positive(d_example, False)
        and is_reduced_positive(d_example, True)
    )
    d_readings = {}
    for name, scale in [("bennett_ghadermarzi_scaled", 1), ("belabas_unscaled", 3)]:
        a, b, c, d = d_example
        b *= scale
        c *= scale
        lo = Fraction((3 * a + b) * c - b * b, 9 * a)
        hi = Fraction((3 * a - b) * c - b * b, 9 * a)
        assert not lo <= d <= hi
        d_readings[name] = dict(printed_lower=str(lo), printed_upper=str(hi), actual_d=d)
    # The printed refinement uses the interval
    # (lower_discriminant, upper_discriminant].
    small = (1, 0, -2, -2)
    assert (
        discriminant(small) == 756
        and is_monic_irreducible(small)
        and is_reduced_positive(small, True)
    )
    a, b, c, d = small
    lower_discriminant, upper_discriminant = 700, 800
    assert lower_discriminant < discriminant(small) <= upper_discriminant
    printed_c_upper = min(
        Fraction(b * b - Fraction(27 * a * a * lower_discriminant, 4), 3 * a), b - 3 * a
    )
    assert not c <= printed_c_upper

    f = (1, 1, 2, 6)
    g = shift_variable(f, 1)
    assert g == (1, 2, 5, 16)
    assert (
        discriminant(f) == discriminant(g) == -216
        and is_monic_irreducible(f)
        and is_monic_irreducible(g)
    )
    assert is_reduced_negative(f, True) and is_reduced_negative(g, True)
    assert not is_reduced_negative(f, False) and not is_reduced_negative(g, False)
    canonical = (1, 0, 1, 2)
    assert shift_variable(canonical, 1) == f and is_reduced_negative(canonical, False)

    h = (1, 2, -1, -2)
    assert discriminant(h) == 2700 and is_monic_irreducible(h)
    assert is_reduced_positive(h, False) and not is_reduced_positive(h, True)

    result = dict(
        d_bound=dict(
            form=d_example,
            discriminant=discriminant(d_example),
            k=-71,
            readings=d_readings,
            point_y_nonnegative=[8, 21],
        ),
        c_refinement=dict(
            form=small,
            discriminant_lower_exclusive=lower_discriminant,
            discriminant_upper_inclusive=upper_discriminant,
            printed_upper=str(printed_c_upper),
            actual_c=c,
        ),
        negative_discriminant_definition=dict(
            equivalent_distinct_forms=[f, g],
            discriminant=-216,
            both_satisfy_printed_definition=True,
            change_of_variables="u -> u + v, v -> v; determinant = 1",
            correct_canonical_form=canonical,
            source_correct_middle="a * d - 9 * b * c",
            printed_middle="3 * (a * d - b * c)",
        ),
        positive_discriminant_boundary=dict(
            form=h,
            discriminant=2700,
            source_reduced=True,
            printed_reduced=False,
            omitted_qualification="Require b < |c| when a = |d|; the printed condition applies it to every P = R case.",
        ),
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
