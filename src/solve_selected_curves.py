#!/usr/bin/env python3
"""Independently enumerate forms and solve their Thue equations with PARI/GP.

Use Belabas’s unscaled-form bounds, converted to (a, 3b, 3c, d), keeping
extra representatives at positive-discriminant reduction boundaries.
Solve the specified curves with PARI/GP, supplied through PATH or --gp.
"""

import argparse
import json
import math
import subprocess
from fractions import Fraction
from functools import reduce
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


def ceil_div(n, d):
    return -((-n) // d)


def ceil_cuberoot(n):
    lo, hi = 0, 1
    while hi**3 < n:
        hi *= 2
    while lo < hi:
        mid = (lo + hi) // 2
        if mid**3 >= n:
            hi = mid
        else:
            lo = mid + 1
    return lo


def enumerate_forms(k):
    D = -108 * k
    T = abs(D)
    A = 0
    denom = 729 if D > 0 else 27
    while (A + 1) ** 4 * denom <= 16 * T:
        A += 1
    result = []
    tested = 0
    for a in range(1, A + 1):
        # Conservative integer upper bound for the square-root expression.
        S = math.isqrt(T if D > 0 else T // 3) + 1
        b_max = (a + 1) // 2 + (math.isqrt(S) + 1 + 2) // 3 + 1
        for b in range(b_max + 1):
            if D > 0:
                c_min = ceil_div(b * b - math.isqrt(T // 81), a)
                c_max = b - a
            else:
                c_min = 1 - b
                W = ceil_cuberoot(ceil_div(T, 4 * a))
                U = Fraction(3 * b * b, a) if a >= 2 * b else Fraction(12 * b - 3 * a, 4)
                c_max = (Fraction(W) + U) // 3
            for c in range(c_min, c_max + 1):
                tested += 1
                # Discriminant equation, divided by -27, is quadratic in d.
                aa = a * a
                bb = 4 * b**3 - 6 * a * b * c
                cc = 4 * a * c**3 - 3 * b * b * c * c - 4 * k
                delta = bb * bb - 4 * aa * cc
                if delta < 0:
                    continue
                rt = math.isqrt(delta)
                if rt * rt != delta:
                    continue
                for num in set([-bb + rt, -bb - rt]):
                    if num % (2 * aa):
                        continue
                    d = num // (2 * aa)
                    if reduce(math.gcd, [a, 3 * b, 3 * c, d]) != 1:
                        continue
                    if D > 0:
                        P = b * b - a * c
                        Q = b * c - a * d
                        R = c * c - b * d
                        if not abs(Q) <= P <= R:
                            continue
                        # Retain extra representatives on the reduction boundary.
                    else:
                        if not d * d - a * a > 3 * (b * d - a * c):
                            continue
                        if (
                            not -((a - 3 * b) ** 2) - 3 * a * c
                            < a * d - 9 * b * c
                            < (a + 3 * b) ** 2 + 3 * a * c
                        ):
                            continue
                        if b == 0 and d <= 0:
                            continue
                    assert (
                        -27
                        * (
                            a * a * d * d
                            - 6 * a * b * c * d
                            - 3 * b * b * c * c
                            + 4 * a * c**3
                            + 4 * b**3 * d
                        )
                        == D
                    )
                    result.append((a, b, c, d))
    return sorted(set(result)), dict(a_max=A, abc_candidates=tested)


def reducible_points(k):
    pts = set()
    four = 4 * abs(k)
    for C in range(-math.isqrt(four), math.isqrt(four) + 1):
        if not C or four % (C * C):
            continue
        q = 4 * C - 4 * k // (C * C)
        if q < 0 or q % 3:
            continue
        B = math.isqrt(q // 3)
        if B * B != q // 3 or B * C % 2:
            continue
        x = B * B - C
        y = abs(B * (3 * C - 2 * B * B) // 2)
        assert y * y - x**3 == k
        pts.add((x, y))
    return pts


def solve_curve(k, folder, gp):
    if k == -432:
        from solve_fermat_exception import solve_exception

        return solve_exception(k)
    fs, bounds = enumerate_forms(k)
    print(f"Found {len(fs):,} candidate forms for k = {k}.", flush=True)
    code = ['print("VERSION ",version());']
    for i, (a, b, c, d) in enumerate(fs):
        P = f"({a})*z^3+({3*b})*z^2+({3*c})*z+({d})"
        code.append(
            f'if(polisirreducible({P}),T=thueinit({P},1);print("F {i} ",thue(T,1)),print("R {i}"));'
        )
    code += ['print("DONE");', "quit;"]
    gpfile = folder / f"curve-{k}.gp"
    gpfile.write_text("\n".join(code) + "\n")
    proc = subprocess.run(
        [gp, "-q", "-s", "128M", str(gpfile)],
        capture_output=True,
        text=True,
        timeout=240,
    )
    (folder / f"curve-{k}.log").write_text(proc.stdout)
    (folder / f"curve-{k}.stderr").write_text(proc.stderr)
    assert proc.returncode == 0 and not proc.stderr.strip(), proc.stderr
    assert proc.stdout.rstrip().endswith("DONE")
    pts = reducible_points(k)
    form_solutions = []
    for line in proc.stdout.splitlines():
        if not line.startswith("F "):
            continue
        _, idx, rest = line.split(" ", 2)
        i = int(idx)
        uvs = json.loads(rest)
        a, b, c, d = fs[i]
        for u, v in uvs:
            assert a * u**3 + 3 * b * u * u * v + 3 * c * u * v * v + d * v**3 == 1
            x = (b * b - a * c) * u * u + (b * c - a * d) * u * v + (c * c - b * d) * v * v
            G = (
                (-a * a * d + 3 * a * b * c - 2 * b**3) * u**3
                + 3 * (-b * b * c - a * b * d + 2 * a * c * c) * u * u * v
                + 3 * (b * c * c - 2 * b * b * d + a * c * d) * u * v * v
                + (-3 * b * c * d + 2 * c**3 + a * d * d) * v**3
            )
            assert G % 2 == 0
            y = abs(G // 2)
            assert y * y - x**3 == k
            pts.add((x, y))
        form_solutions.append(dict(form=[a, b, c, d], solutions=uvs))
    return dict(
        k=k,
        enumeration_bounds=bounds,
        candidate_forms=fs,
        irreducible_forms_solved=len(form_solutions),
        thue_solutions=form_solutions,
        points_y_nonnegative=sorted(pts),
        affine_point_count=sum(1 if y == 0 else 2 for x, y in pts),
        thueinit_flag=1,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("coefficients", nargs="+", type=int)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--gp", default="gp")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for k in args.coefficients:
        result = solve_curve(k, args.output_dir, args.gp)
        (args.output_dir / f"curve-{k}.json").write_text(json.dumps(result, indent=2) + "\n")
        print(
            json.dumps(
                {
                    key: v
                    for key, v in result.items()
                    if key not in ["thue_solutions", "candidate_forms"]
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
