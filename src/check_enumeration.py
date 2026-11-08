#!/usr/bin/env python3
"""Compare interval form sets with independent per-coefficient enumeration."""

import argparse
import concurrent.futures
import json
import random
import time
from pathlib import Path

from common import sha256_file
from solve_selected_curves import enumerate_forms

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


def check(item):
    k, expected = item
    actual, _ = enumerate_forms(k)
    assert set(actual) == set(expected), (k, actual, expected)
    return dict(k=k, forms=len(actual))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enumeration-dir", type=Path, required=True)
    ap.add_argument("--sample-size", type=int, default=3000)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seed", type=int, default=20261103)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    start = time.monotonic()
    meta = json.loads((args.enumeration_dir / "enumeration.json").read_text())
    lo, hi = meta["lower_inclusive"], meta["upper_inclusive"]
    candidates = [k for k in range(lo, hi + 1) if k]
    ks = set(random.Random(args.seed).sample(candidates, min(args.sample_size, len(candidates))))
    # Include the coefficient with the large omitted point used in the note.
    ks.update(k for k in [lo, lo + 1, hi - 1, hi, 1436108] if lo <= k <= hi and k)
    expected = {k: [] for k in ks}
    with (args.enumeration_dir / "forms.tsv").open() as f:
        for line in f:
            k, *form = map(int, line.split())
            if k in expected:
                expected[k].append(tuple(form))
    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for r in pool.map(check, sorted(expected.items()), chunksize=5):
            results.append(r)
            if len(results) % 100 == 0:
                print(
                    f"Checked {len(results):,} coefficients in {time.monotonic() - start:.2f} seconds.",
                    flush=True,
                )

    result = dict(
        all_identical=True,
        seed=args.seed,
        random_sample_size=args.sample_size,
        coefficients_checked=len(results),
        elapsed_seconds=time.monotonic() - start,
        forms_sha256=sha256_file(args.enumeration_dir / "forms.tsv"),
        checks=results,
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Form sets agree for all {len(results):,} coefficients.", flush=True)


if __name__ == "__main__":
    main()
