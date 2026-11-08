#!/usr/bin/env python3
"""Independently compare the full interval export, published data and box search."""

import argparse
import csv
import json
from pathlib import Path

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result-dir", type=Path, required=True)
    ap.add_argument("--enumeration-dir", type=Path, required=True)
    ap.add_argument("--published-rows", type=Path, required=True)
    ap.add_argument("--known-missing", type=Path, required=True)
    ap.add_argument(
        "--output",
        type=Path,
        help="Write the report here instead of result-dir/comparison.json.",
    )
    args = ap.parse_args()
    root = args.result_dir
    ed = args.enumeration_dir
    summary = json.loads((root / "summary.json").read_text())
    assert summary["all_complete"] and not summary["errors"]
    lo, hi = summary["lower_inclusive"], summary["upper_inclusive"]

    def read_csv(p):
        with p.open() as f:
            reader = csv.DictReader(f)
            assert reader.fieldnames == ["k", "x", "y"]
            rows = [(int(r["k"]), int(r["x"]), int(r["y"])) for r in reader]
        assert len(rows) == len(set(rows)), f"Duplicate point rows in {p}."
        return set(rows)

    computed = read_csv(root / "points.csv")
    assert len(computed) == summary["affine_points"]
    assert all(lo <= k <= hi and y * y - x * x * x == k for k, x, y in computed)
    published = set()
    present = set()
    with args.published_rows.open() as f:
        for line in f:
            k, ps = json.loads(line)
            if lo <= k <= hi:
                present.add(k)
                published.update((k, x, y) for x, y in ps)
    assert not published - computed, "Some published points were not recovered."
    missing = computed - published
    assert missing == read_csv(root / "missing-points.csv")
    recorded = {
        (r["k"], x, yy)
        for r in summary["published_differences"]
        for x, y in r["computed_only"]
        for yy in ([0] if y == 0 else [-y, y])
    }
    assert recorded == missing and not any(
        r["published_only"] for r in summary["published_differences"]
    )
    known = {p for p in read_csv(args.known_missing) if lo <= p[0] <= hi}
    assert known <= missing, "A previously established omission was not reproduced."
    boxmeta = json.loads((ed / "box.json").read_text())
    assert (boxmeta["lower_inclusive"], boxmeta["upper_inclusive"]) == (lo, hi)
    box = set()
    for line in (ed / "box.tsv").open():
        k, x, y = map(int, line.split())
        assert y >= 0 and y * y - x * x * x == k
        assert (k, x, y) not in box
        box.add((k, x, y))
    assert len(box) == boxmeta["points_y_nonnegative"]
    computed_box = {(k, x, y) for k, x, y in computed if y >= 0 and x <= boxmeta["x_max"]}
    assert box == computed_box, "The independent coordinate-box search disagrees."
    affected = {k for k, x, y in missing}
    result = dict(
        coefficients=summary["coefficients"],
        affine_points=len(computed),
        published_affine_points=len(published),
        published_points_not_recovered=0,
        missing_affine_points=len(missing),
        missing_points_y_nonnegative=sum(y >= 0 for k, x, y in missing),
        affected_coefficients=len(affected),
        absent_published_rows=len(affected - present),
        incomplete_published_rows=len(affected & present),
        earlier_missing_affine_points_reproduced=len(known),
        additional_missing_affine_points=len(missing - known),
        additional_missing_points=sorted(missing - known),
        independent_box=dict(**boxmeta, exact_agreement=True),
    )
    (args.output or root / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
