"""Cross-check the complete tables, published snapshots and derived data."""

import csv
import gzip
import json
from collections import Counter
from itertools import chain


def require(condition, message):
    if not condition:
        raise ValueError(message)


def published_rows(root):
    # Negative source rows run from -1 downwards; reverse their text records.
    # Buffer this smaller snapshot and stream the full table.
    previous = -10000001
    for sign in ["negative", "positive"]:
        with gzip.open(root / f"data/published/{sign}-rows.jsonl.gz", "rt") as stream:
            lines = reversed(list(stream)) if sign == "negative" else stream
            for line in lines:
                k, points = json.loads(line)
                require(
                    previous < k and 0 < abs(k) <= 10000000,
                    "Published rows are out of order or outside the coefficient range.",
                )
                require(
                    (k < 0) == (sign == "negative"),
                    "Published coefficient has the wrong sign.",
                )
                previous = k
                signed = set(map(tuple, points))
                require(
                    len(signed) == len(points) > 0,
                    f"Published row is empty or contains duplicate points at k = {k}.",
                )
                yield k, signed


def verify_data(root):
    stats = json.loads((root / "data/statistics.json").read_text())
    corrections = json.loads((root / "data/corrected-curves.json").read_text())
    corrected = {r["k"]: r for r in corrections}
    require(
        len(corrected) == len(corrections),
        "Corrected entries contain a duplicate coefficient.",
    )
    with (root / "data/missing-points.csv").open() as stream:
        reader = csv.reader(stream)
        require(
            next(reader) == ["k", "x", "y"],
            "Unexpected columns in the missing-point CSV.",
        )
        missing_rows = [tuple(map(int, row)) for row in reader]
    expected_missing = set(missing_rows)
    require(
        len(expected_missing) == len(missing_rows),
        "The missing-point CSV contains a duplicate point.",
    )
    found_missing = set()
    affected = set()
    absent = 0
    counts = {s: Counter() for s in ["negative", "positive"]}
    published = iter(published_rows(root))
    baseline = next(published, None)
    with (
        gzip.open(root / "data/curves.jsonl.gz", "rt") as curves,
        gzip.open(root / "data/points.csv.gz", "rt") as points,
    ):
        reader = csv.reader(points)
        require(next(reader) == ["k", "x", "y"], "Unexpected columns in the full-point CSV.")
        for expected_k in chain(range(-10000000, 0), range(1, 10000001)):
            line = next(curves, None)
            require(line is not None, f"Missing curve record at k = {expected_k}.")
            row = json.loads(line)
            k = row["k"]
            require(
                k == expected_k and row["status"] == "complete",
                f"Curve order or completion status is invalid at k = {expected_k}.",
            )
            normalized = row["points_y_nonnegative"]
            require(
                len(normalized) == len(set(map(tuple, normalized))),
                f"Curve record contains a duplicate point at k = {k}.",
            )
            require(
                all(
                    type(x) is int and type(y) is int and y >= 0 and y * y == x * x * x + k
                    for x, y in normalized
                ),
                f"Curve record contains an invalid point at k = {k}.",
            )
            signed = {(x, yy) for x, y in normalized for yy in {-y, y}}
            require(
                len(signed) == row["affine_point_count"],
                f"Curve point count is incorrect at k = {k}.",
            )
            for x, y in sorted(signed):
                record = next(reader, None)
                require(
                    record is not None and tuple(map(int, record)) == (k, x, y),
                    f"CSV and JSONL points disagree at k = {k}.",
                )
            original = set()
            require(
                baseline is None or baseline[0] >= k,
                f"A published row was skipped before k = {k}.",
            )
            if baseline is not None and baseline[0] == k:
                original = baseline[1]
                baseline = next(published, None)
            require(original <= signed, f"A published point was not recovered at k = {k}.")
            extra = signed - original
            if extra:
                affected.add(k)
                absent += not original
                found_missing.update((k, x, y) for x, y in extra)
                correction = corrected.get(k)
                require(correction is not None, f"A correction is missing at k = {k}.")
                require(
                    correction["published_affine_points"] == len(original)
                    and correction["corrected_affine_points"] == len(signed)
                    and correction["corrected_points_y_nonnegative"] == normalized,
                    f"The corrected entry disagrees at k = {k}.",
                )
            counts["negative" if k < 0 else "positive"][len(signed)] += 1
        require(
            next(curves, None) is None and next(reader, None) is None,
            "The full tables contain extra records.",
        )
    require(baseline is None, "Published rows remain unmatched.")
    require(
        found_missing == expected_missing and affected == set(corrected),
        "The correction files disagree with the full comparison.",
    )
    require(
        (len(found_missing), len(affected), absent) == (783, 389, 244),
        "The omission totals disagree.",
    )
    for sign, histogram in counts.items():
        actual = dict(
            curves=sum(histogram.values()),
            nonempty_curves=sum(histogram.values()) - histogram[0],
            affine_points=sum(n * c for n, c in histogram.items()),
            histogram={str(n): c for n, c in sorted(histogram.items())},
            empty_curves=histogram[0],
        )
        require(
            actual == stats["by_sign"][sign],
            f"Statistics disagree for {sign} coefficients.",
        )
    total = sum(sum(n * c for n, c in h.items()) for h in counts.values())
    require(
        stats["scope"] == "0 < |k| <= 10000000",
        "The statistics have an unexpected scope.",
    )
    require(
        stats["coefficients"] == 20000000 and stats["affine_points"] == total,
        "The aggregate totals disagree.",
    )
    print(
        f"Checked all 20,000,000 curves and {total:,} signed points. "
        "Found 783 omissions on 389 curves. All statistics agree."
    )
