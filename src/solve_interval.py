#!/usr/bin/env python3
"""Solve every enumerated form in batches and export points for each coefficient.

Read the TSV files from enumerate_forms.cpp. Record each form as F (solved),
R (reducible), A (proved algebraically), or E (failed). Retain errors and mark
affected curves as incomplete.
"""

import argparse
import collections
import concurrent.futures
import csv
import hashlib
import json
import subprocess
import time
from pathlib import Path

from common import sha256_file
from solve_fermat_exception import solve_exception

if not __debug__:
    raise RuntimeError("Run without -O: Verification requires assertions.")


# Used only to prevent incompatible checkpoint reuse during --resume.
SOLVER_SOURCES = (
    "solve_interval.py",
    "enumerate_forms.cpp",
    "solve_fermat_exception.py",
    "solve_selected_curves.py",
    "common.py",
)


def point(f, uv, k):
    a, b, c, d = f
    u, v = uv
    assert a * u**3 + 3 * b * u * u * v + 3 * c * u * v * v + d * v**3 == 1
    x = (b * b - a * c) * u * u + (b * c - a * d) * u * v + (c * c - b * d) * v * v
    g = (
        (-a * a * d + 3 * a * b * c - 2 * b**3) * u**3
        + 3 * (-b * b * c - a * b * d + 2 * a * c * c) * u * u * v
        + 3 * (b * c * c - 2 * b * b * d + a * c * d) * u * v * v
        + (-3 * b * c * d + 2 * c**3 + a * d * d) * v**3
    )
    assert g % 2 == 0
    y = abs(g // 2)
    assert y * y - x * x * x == k
    return [x, y]


def solve_job(item):
    job_id, rows, directory, gp, resume = item
    folder = Path(directory) / "jobs"
    folder.mkdir(exist_ok=True)
    stem = folder / f"job-{job_id:06d}"
    code = ['print("VERSION ",version());']
    by_id = {i: (k, f) for i, k, f in rows}
    assert len(by_id) == len(rows)
    for i, k, (a, b, c, d) in rows:
        if k == -432:
            assert [a, b, c, d] == [1, 3, -3, -9]
            code.append(f'print("A {i} [[1,0]]");')
        else:
            p = f"({a})*z^3+({3*b})*z^2+({3*c})*z+({d})"
            code.append(
                f'setrand(1);P={p};if(polisirreducible(P),iferr(S=thue(thueinit(P,1),1);print("F {i} ",S),err,print("E {i} ",err)),print("R {i}"));'
            )
    code += ['print("DONE");', "quit;"]
    program = "\n".join(code) + "\n"
    ph = hashlib.sha256(program.encode()).hexdigest()
    checkpoint = stem.with_suffix(".json")
    if resume and checkpoint.exists():
        result = json.loads(checkpoint.read_text())
        assert result["program_sha256"] == ph and result["log_sha256"] == sha256_file(
            stem.with_suffix(".log")
        )
        assert result["stderr_sha256"] == sha256_file(stem.with_suffix(".stderr"))
        return result
    stem.with_suffix(".gp").write_text(program)
    started = time.perf_counter()
    proc = subprocess.run(
        [gp, "-q", "-s", "128M", str(stem.with_suffix(".gp"))],
        capture_output=True,
        text=True,
        timeout=3600,
    )
    stem.with_suffix(".log").write_text(proc.stdout)
    stem.with_suffix(".stderr").write_text(proc.stderr)
    assert proc.returncode == 0 and not proc.stderr.strip(), f"Job {job_id}: {proc.stderr}"
    lines = proc.stdout.splitlines()
    assert lines[0].startswith("VERSION ") and lines[-1] == "DONE"
    records = []
    seen = set()
    for line in lines[1:-1]:
        parts = line.split(" ", 2)
        kind = parts[0]
        i = int(parts[1])
        assert i in by_id and i not in seen
        seen.add(i)
        k, f = by_id[i]
        record = dict(i=i, k=k, kind=kind)
        if kind in ["F", "A"]:
            solutions = json.loads(parts[2])
            assert len(solutions) == len(set(map(tuple, solutions)))
            record.update(solutions=solutions, points=[point(f, uv, k) for uv in solutions])
            if kind == "A":
                assert k == -432 and solutions == [[1, 0]]
        elif kind == "E":
            record["error"] = parts[2]
        else:
            assert kind == "R" and len(parts) == 2
        records.append(record)
    assert seen == set(by_id)
    result = dict(
        job=job_id,
        seconds=time.perf_counter() - started,
        program_sha256=ph,
        log_sha256=sha256_file(stem.with_suffix(".log")),
        stderr_sha256=sha256_file(stem.with_suffix(".stderr")),
        records=records,
    )
    temporary = checkpoint.with_suffix(".tmp")
    temporary.write_text(json.dumps(result, separators=(",", ":")) + "\n")
    temporary.replace(checkpoint)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--enumeration-dir", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--chunk-size", type=int, default=1000)
    ap.add_argument("--gp", default="gp")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument(
        "--compare-curves",
        type=Path,
        help="Compare with an independently validated small-range JSONL file.",
    )
    ap.add_argument("--published-rows", type=Path, nargs="*")
    args = ap.parse_args()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    assert 1 <= args.workers <= 32 and 1 <= args.chunk_size <= 10000
    enumeration = json.loads((args.enumeration_dir / "enumeration.json").read_text())
    lo, hi = enumeration["lower_inclusive"], enumeration["upper_inclusive"]
    if lo <= -432 <= hi:
        exception = solve_exception(-432)
        assert exception["thue_solutions"][0]["solutions"] == [(1, 0)]
    config = dict(
        enumeration=enumeration,
        forms_sha256=sha256_file(args.enumeration_dir / "forms.tsv"),
        reducible_sha256=sha256_file(args.enumeration_dir / "reducible.tsv"),
        chunk_size=args.chunk_size,
        source_hashes={
            name: sha256_file(Path(__file__).with_name(name)) for name in SOLVER_SOURCES
        },
        gp_version=subprocess.check_output([args.gp, "--version-short"], text=True).strip(),
    )
    cp = out / "run-config.json"
    if cp.exists():
        assert args.resume and json.loads(cp.read_text()) == config
    else:
        cp.write_text(json.dumps(config, indent=2) + "\n")
    points = collections.defaultdict(set)
    form_counts = collections.Counter()
    irred_counts = collections.Counter()
    errors = []
    kinds = collections.Counter()
    jobs = 0
    processed = 0
    for line in (args.enumeration_dir / "reducible.tsv").open():
        k, x, y, B, C = map(int, line.split())
        assert lo <= k <= hi and k and y >= 0
        assert 4 * k == C * C * (4 * C - 3 * B * B) and B * C % 2 == 0
        assert x == B * B - C and y == abs(B * (3 * C - 2 * B * B) // 2) and y * y - x * x * x == k
        points[k].add((x, y))
    started = time.perf_counter()

    def chunks():
        rows = []
        job = 0
        for i, line in enumerate((args.enumeration_dir / "forms.tsv").open()):
            k, a, b, c, d = map(int, line.split())
            assert lo <= k <= hi and k
            assert (
                a * a * d * d
                - 6 * a * b * c * d
                - 3 * b * b * c * c
                + 4 * a * c * c * c
                + 4 * b * b * b * d
                == 4 * k
            )
            form_counts[k] += 1
            rows.append((i, k, [a, b, c, d]))
            if len(rows) == args.chunk_size:
                yield job, rows, str(out), args.gp, args.resume
                job += 1
                rows = []
        if rows:
            yield job, rows, str(out), args.gp, args.resume

    iterator = iter(chunks())
    with (out / "progress.jsonl").open("a" if args.resume else "w") as progress:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
            pending = set()

            def submit_one():
                try:
                    pending.add(pool.submit(solve_job, next(iterator)))
                    return True
                except StopIteration:
                    return False

            for _ in range(2 * args.workers):
                submit_one()
            while pending:
                done, pending = concurrent.futures.wait(
                    pending, return_when=concurrent.futures.FIRST_COMPLETED
                )
                for future in done:
                    result = future.result()
                    jobs += 1
                    for r in result["records"]:
                        processed += 1
                        kinds[r["kind"]] += 1
                        if r["kind"] in ["F", "A"]:
                            irred_counts[r["k"]] += 1
                            points[r["k"]].update(map(tuple, r["points"]))
                        elif r["kind"] == "E":
                            errors.append(r)
                    entry = dict(
                        jobs=jobs,
                        forms=processed,
                        total_forms=enumeration["candidate_forms"],
                        errors=len(errors),
                        elapsed_seconds=round(time.perf_counter() - started, 2),
                    )
                    progress.write(json.dumps(entry) + "\n")
                    progress.flush()
                    if jobs % 10 == 0 or len(errors):
                        print(json.dumps(entry), flush=True)
                    submit_one()
    assert processed == enumeration["candidate_forms"] == sum(form_counts.values())
    compute_seconds = time.perf_counter() - started
    published = {}
    for path in args.published_rows or []:
        for line in path.open():
            k = int(line[1 : line.index(",")])
            if lo <= k <= hi:
                key, ps = json.loads(line)
                assert key == k
                assert all(y * y - x * x * x == k for x, y in ps)
                published[k] = {(x, abs(y)) for x, y in ps}
    expected = {}
    if args.compare_curves:
        for line in args.compare_curves.open():
            r = json.loads(line)
            expected[r["k"]] = set(map(tuple, r["points_y_nonnegative"]))
        assert set(expected) == {k for k in range(lo, hi + 1) if k}
    bad = {r["k"] for r in errors}
    mismatches = []
    differences = []
    signed = normalized = empty = 0
    total = 0
    with (
        (out / "curves.jsonl").open("w") as curves,
        (out / "points.csv").open("w", newline="") as pf,
        (out / "missing-points.csv").open("w", newline="") as mf,
    ):
        pw = csv.writer(pf)
        mw = csv.writer(mf)
        pw.writerow(["k", "x", "y"])
        mw.writerow(["k", "x", "y"])
        for k in range(lo, hi + 1):
            if not k:
                continue
            total += 1
            ps = points.get(k, set())
            n = sum(1 if y == 0 else 2 for x, y in ps)
            curves.write(
                json.dumps(
                    dict(
                        k=k,
                        status="incomplete" if k in bad else "complete",
                        points_y_nonnegative=sorted(ps),
                        affine_point_count=n,
                        candidate_forms=form_counts[k],
                        irreducible_forms_solved=irred_counts[k],
                    ),
                    separators=(",", ":"),
                )
                + "\n"
            )
            for x, y in sorted(ps):
                for yy in ([0] if not y else [-y, y]):
                    pw.writerow([k, x, yy])
            signed += n
            normalized += len(ps)
            empty += not ps
            if expected and ps != expected[k]:
                mismatches.append(
                    dict(
                        k=k,
                        computed_only=sorted(ps - expected[k]),
                        expected_only=sorted(expected[k] - ps),
                    )
                )
            if args.published_rows and ps != published.get(k, set()):
                old = published.get(k, set())
                new = ps - old
                differences.append(
                    dict(
                        k=k,
                        computed_only=sorted(new),
                        published_only=sorted(old - ps),
                        published_row_absent=k not in published,
                    )
                )
                for x, y in sorted(new):
                    for yy in ([0] if not y else [-y, y]):
                        mw.writerow([k, x, yy])
    summary = dict(
        lower_inclusive=lo,
        upper_inclusive=hi,
        coefficients=total,
        completed_coefficients=total - len(bad),
        jobs=jobs,
        candidate_forms=processed,
        form_status_counts=dict(kinds),
        errors=errors,
        comparison_mismatches=mismatches,
        published_differences=differences,
        affine_points=signed,
        points_y_nonnegative=normalized,
        zero_solution_curves=empty,
        solving_wall_seconds=compute_seconds,
        all_complete=not errors,
        all_comparisons_match=not mismatches,
    )
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k != "published_differences"},
            indent=2,
        ),
        flush=True,
    )
    if errors or mismatches:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
