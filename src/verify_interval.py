#!/usr/bin/env python3
"""Check raw solver results by exact arithmetic and verify the interval exports.

Check that the jobs cover every enumerated form and reconstruct the exported
points. Completeness relies on the reduction bounds and PARI’s Thue solver.
"""

import argparse
import collections
import csv
import json
import math
from pathlib import Path

from common import sha256_file

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


def image(f, uv, k):
    a, b, c, d = f
    u, v = uv
    assert a * u**3 + 3 * b * u * u * v + 3 * c * u * v * v + d * v**3 == 1
    # Recover the point independently using derivatives of the form and Hessian.
    P = b * b - a * c
    Q = b * c - a * d
    R = c * c - b * d
    x = P * u * u + Q * u * v + R * v * v
    Fu = 3 * a * u * u + 6 * b * u * v + 3 * c * v * v
    Fv = 3 * b * u * u + 6 * c * u * v + 3 * d * v * v
    Hu = 2 * P * u + Q * v
    Hv = Q * u + 2 * R * v
    # The normalized Jacobian is (Fv * Hu - Fu * Hv) / 3.
    numerator = Fv * Hu - Fu * Hv
    assert numerator % 6 == 0
    y = abs(numerator // 6)
    assert y * y - x * x * x == k
    return x, y


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enumeration-dir", type=Path, required=True)
    ap.add_argument("--result-dir", type=Path, required=True)
    ap.add_argument("--small-reference", type=Path)
    ap.add_argument("--affected-reference", type=Path)
    ap.add_argument(
        "--reducibility-check",
        type=Path,
        help="Check against a complete, independent rational-root classification JSON file.",
    )
    args = ap.parse_args()
    ed = args.enumeration_dir
    root = args.result_dir
    config = json.loads((root / "run-config.json").read_text())
    summary = json.loads((root / "summary.json").read_text())
    assert (
        summary["all_complete"] and not summary["errors"] and not summary["comparison_mismatches"]
    )
    assert (
        sha256_file(ed / "forms.tsv") == config["forms_sha256"]
        and sha256_file(ed / "reducible.tsv") == config["reducible_sha256"]
    )
    independent_reducible = None
    if args.reducibility_check:
        classification = json.loads(args.reducibility_check.read_text())
        assert (
            classification["forms_sha256"] == config["forms_sha256"]
            and classification["forms_checked"] == summary["candidate_forms"]
        )
        independent_reducible = set(classification["reducible_indices"])
    lo, hi = summary["lower_inclusive"], summary["upper_inclusive"]
    pts = collections.defaultdict(set)
    counts = collections.Counter()
    irred = collections.Counter()
    nr = 0
    for line in (ed / "reducible.tsv").open():
        k, x, y, B, C = map(int, line.split())
        assert lo <= k <= hi and k and y >= 0
        assert 4 * k == C * C * (4 * C - 3 * B * B) and B * C % 2 == 0
        assert x == B * B - C and y == abs(B * (3 * C - 2 * B * B) // 2) and y * y - x * x * x == k
        pts[k].add((x, y))
        nr += 1
    assert nr == config["enumeration"]["reducible_records"]
    form_stream = (ed / "forms.tsv").open()
    next_id = 0
    kinds = collections.Counter()
    uvcount = 0
    jobs = sorted((root / "jobs").glob("job-*.json"))
    for job_index, jpath in enumerate(jobs):
        assert jpath.stem == f"job-{job_index:06d}"
        record = json.loads(jpath.read_text())
        assert record["job"] == job_index
        assert sha256_file(jpath.with_suffix(".gp")) == record["program_sha256"]
        assert sha256_file(jpath.with_suffix(".log")) == record["log_sha256"]
        assert sha256_file(jpath.with_suffix(".stderr")) == record["stderr_sha256"]
        assert not jpath.with_suffix(".stderr").read_text().strip()
        lines = jpath.with_suffix(".log").read_text().splitlines()
        assert lines[0].startswith("VERSION ") and lines[-1] == "DONE"
        assert len(lines) - 2 == len(record["records"])
        program = ['print("VERSION ",version());']
        for line, saved in zip(lines[1:-1], record["records"]):
            raw = form_stream.readline()
            assert raw
            k, a, b, c, d = map(int, raw.split())
            f = (a, b, c, d)
            assert lo <= k <= hi and k and math.gcd(a, 3 * b, 3 * c, d) == 1
            assert (
                a * a * d * d
                - 6 * a * b * c * d
                - 3 * b * b * c * c
                + 4 * a * c * c * c
                + 4 * b * b * b * d
                == 4 * k
            )
            counts[k] += 1
            parts = line.split(" ", 2)
            kind = parts[0]
            i = int(parts[1])
            assert i == next_id == saved["i"]
            assert saved["kind"] == kind and saved["k"] == k
            if independent_reducible is not None:
                assert (kind == "R") == (i in independent_reducible)
            next_id += 1
            kinds[kind] += 1
            if k == -432:
                assert f == (1, 3, -3, -9) and kind == "A"
                program.append(f'print("A {i} [[1,0]]");')
            else:
                p = f"({a})*z^3+({3*b})*z^2+({3*c})*z+({d})"
                program.append(
                    f'setrand(1);P={p};if(polisirreducible(P),iferr(S=thue(thueinit(P,1),1);print("F {i} ",S),err,print("E {i} ",err)),print("R {i}"));'
                )
            if kind in ["F", "A"]:
                uvs = json.loads(parts[2])
                assert uvs == saved["solutions"] and len(uvs) == len(set(map(tuple, uvs)))
                images = [image(f, uv, k) for uv in uvs]
                assert images == list(map(tuple, saved["points"]))
                if kind == "A":
                    assert uvs == [[1, 0]] and images == [(12, 36)]
                pts[k].update(images)
                irred[k] += 1
                uvcount += len(uvs)
            else:
                assert kind == "R" and len(parts) == 2
        program += ['print("DONE");', "quit;"]
        assert jpath.with_suffix(".gp").read_text() == "\n".join(program) + "\n"
        if (job_index + 1) % 100 == 0:
            print(
                f"Verified {job_index + 1:,} jobs covering {next_id:,} forms.",
                flush=True,
            )
    assert not form_stream.readline()
    form_stream.close()
    assert next_id == summary["candidate_forms"] and dict(kinds) == summary["form_status_counts"]
    small = {}
    if args.small_reference:
        for line in args.small_reference.open():
            r = json.loads(line)
            small[r["k"]] = r["points_y_nonnegative"]
    affected = {}
    if args.affected_reference:
        for r in json.loads(args.affected_reference.read_text()):
            if lo <= r["k"] <= hi:
                affected[r["k"]] = r["corrected_points_y_nonnegative"]
    total = signed = normalized = empty = affected_checked = 0
    expected_k = lo
    with (root / "curves.jsonl").open() as curves:
        for line in curves:
            if expected_k == 0:
                expected_k += 1
            r = json.loads(line)
            k = r["k"]
            assert k == expected_k and r["status"] == "complete"
            expected_k += 1
            assert r["points_y_nonnegative"] == [list(p) for p in sorted(pts.get(k, set()))]
            assert r["candidate_forms"] == counts[k] and r["irreducible_forms_solved"] == irred[k]
            n = sum(1 if y == 0 else 2 for x, y in pts.get(k, set()))
            assert r["affine_point_count"] == n
            if small:
                assert r["points_y_nonnegative"] == small[k]
            if k in affected:
                assert r["points_y_nonnegative"] == affected[k]
                affected_checked += 1
            total += 1
            signed += n
            normalized += len(pts.get(k, set()))
            empty += not pts.get(k, set())
    assert expected_k == hi + 1 and total == summary["coefficients"]
    assert (signed, normalized, empty) == (
        summary["affine_points"],
        summary["points_y_nonnegative"],
        summary["zero_solution_curves"],
    )
    csv_points = set()
    with (root / "points.csv").open() as stream:
        reader = csv.DictReader(stream)
        assert reader.fieldnames == ["k", "x", "y"]
        for r in reader:
            k, x, y = map(int, [r["k"], r["x"], r["y"]])
            assert y * y - x * x * x == k
            assert (x, abs(y)) in pts[k] and (k, x, y) not in csv_points
            csv_points.add((k, x, y))
    assert len(csv_points) == signed
    # Sign symmetry and cardinality ensure that the CSV contains every point.
    assert all((k, x, -y) in csv_points for k, x, y in csv_points)
    result = dict(
        coefficients=total,
        affine_points=signed,
        candidate_forms=next_id,
        raw_jobs=len(jobs),
        form_status_counts=dict(kinds),
        thue_solution_pairs_checked=uvcount,
        reducible_records=nr,
        earlier_affected_curves_reproduced=affected_checked,
        small_reference_compared=bool(small),
        independent_reducibility_checked=independent_reducible is not None,
        verification_method="Recover points from the raw GP output using the Jacobian in derivative form. Check job coverage, indices, discriminants, equations, curve records, and the CSV export.",
    )
    (root / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
