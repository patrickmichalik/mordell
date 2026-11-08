#!/usr/bin/env python3
"""Count curves and signed points in the complete, ordered Mordell table."""

import argparse
import collections
import json
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("curves", type=Path)
    ap.add_argument("output", type=Path)
    args = ap.parse_args()
    counts = {
        s: dict(
            curves=0,
            nonempty_curves=0,
            affine_points=0,
            histogram=collections.Counter(),
        )
        for s in ["positive", "negative"]
    }
    previous = -10000001
    with args.curves.open("rb") as f:
        for line in f:
            r = json.loads(line)
            k = r["k"]
            assert k == (previous + 1 if previous != -1 else 1) and r["status"] == "complete"
            previous = k
            count = sum(1 if y == 0 else 2 for x, y in r["points_y_nonnegative"])
            assert count == r["affine_point_count"]
            c = counts["positive" if k > 0 else "negative"]
            c["curves"] += 1
            c["nonempty_curves"] += bool(count)
            c["affine_points"] += count
            c["histogram"][count] += 1
    assert previous == 10000000
    for c in counts.values():
        c["empty_curves"] = c["curves"] - c["nonempty_curves"]
        c["histogram"] = {str(k): v for k, v in sorted(c["histogram"].items())}
    report = dict(
        scope="0 < |k| <= 10000000",
        coefficients=sum(c["curves"] for c in counts.values()),
        affine_points=sum(c["affine_points"] for c in counts.values()),
        by_sign=counts,
    )
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
