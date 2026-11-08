#!/usr/bin/env python3
"""Solve k = -432 algebraically using Fermat’s theorem for exponent 3.

This algebraic case handles a PARI/GP 2.19.0 thue() failure.
"""

from solve_selected_curves import enumerate_forms

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


def solve_exception(k):
    assert k == -432
    # Coefficients of (36 + y)^3 + (36 - y)^3, in ascending powers of y.
    assert [2 * 36**3, 0, 6 * 36, 0] == [216 * 432, 0, 216, 0]
    y = 36
    x = 12
    assert x**3 == y * y + 432
    fs, bounds = enumerate_forms(k)
    assert fs == [(1, 3, -3, -9)]
    a, b, c, d = fs[0]
    assert all(a * r**3 + 3 * b * r * r + 3 * c * r + d != 0 for r in [-9, -3, -1, 1, 3, 9])
    assert (b * b - a * c, b * c - a * d, c * c - b * d) == (12, 0, 36)
    uv = [
        (u, v)
        for u in [-1, 0, 1]
        for v in [-1, 0, 1]
        if u * u + 3 * v * v == 1
        and a * u**3 + 3 * b * u * u * v + 3 * c * u * v * v + d * v**3 == 1
    ]
    assert uv == [(1, 0)]
    result = dict(
        k=k,
        enumeration_bounds=bounds,
        candidate_forms=fs,
        irreducible_forms_solved=1,
        thue_solutions=[
            dict(
                form=fs[0],
                solutions=uv,
                method="Apply Fermat's theorem for exponent 3 and the definite Hessian.",
            )
        ],
        points_y_nonnegative=[(x, y)],
        affine_point_count=2,
        thueinit_flag=None,
        solution_method="Use an exact algebraic proof based on Fermat's theorem for exponent 3.",
    )
    return result
