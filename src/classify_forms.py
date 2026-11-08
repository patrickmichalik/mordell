#!/usr/bin/env python3
"""Exact cubic reducibility test over Q, independent of PARI.

A cubic is reducible over Q if and only if it has a rational root. For leading
coefficient a and nonzero constant d, a reduced root m / n satisfies m | d
and n | a. A rootless reduction modulo a prime not dividing a proves
irreducibility and is an optional shortcut; otherwise, every rational-root
candidate is checked exactly.
"""

import argparse
import functools
import json
import math
import time
from pathlib import Path

from common import sha256_file

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


@functools.lru_cache(maxsize=None)
def divisors(n):
    assert n > 0
    small = [d for d in range(1, math.isqrt(n) + 1) if n % d == 0]
    return tuple(sorted(set(small + [n // d for d in small])))


@functools.lru_cache(maxsize=None)
def no_mod_root(a, b, c, d, p):
    return all((a * r * r * r + b * r * r + c * r + d) % p for r in range(p))


def rational_root(f):
    a, b, c, d = f
    b *= 3
    c *= 3
    if d == 0:
        return (0, 1)
    for p in [2, 3, 5, 7, 11]:
        if a % p and no_mod_root(a % p, b % p, c % p, d % p, p):
            return None
    for n in divisors(abs(a)):
        for m in divisors(abs(d)):
            if math.gcd(m, n) != 1:
                continue
            for numerator in [m, -m]:
                if (
                    a * numerator**3
                    + b * numerator * numerator * n
                    + c * numerator * n * n
                    + d * n**3
                    == 0
                ):
                    return numerator, n
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--forms", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    started = time.monotonic()
    reducible = []
    roots = []
    count = 0
    with args.forms.open() as stream:
        for i, line in enumerate(stream):
            k, *f = map(int, line.split())
            root = rational_root(tuple(f))
            count += 1
            if root is not None:
                reducible.append(i)
                roots.append(root)
            if count % 100000 == 0:
                print(
                    f"Classified {count:,} forms; {len(reducible):,} are reducible.",
                    flush=True,
                )
    result = dict(
        forms_checked=count,
        reducible_forms=len(reducible),
        irreducible_forms=count - len(reducible),
        reducible_indices=reducible,
        rational_roots=roots,
        forms_sha256=sha256_file(args.forms),
        elapsed_seconds=time.monotonic() - started,
        method="Apply the rational-root theorem with rootless finite-field shortcuts.",
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"Classification complete: {len(reducible):,} of {count:,} forms are reducible.",
        flush=True,
    )


if __name__ == "__main__":
    main()
