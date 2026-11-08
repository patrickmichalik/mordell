# Corrected integral-point data for Mordell curves

Data and code accompanying the
[research note](https://patrickmichalik.com/notes/corrected-integral-point-data-for-mordell-curves).
Please cite the note when using these data in a publication. See [`LICENSE.md`](LICENSE.md) for
reuse terms.

## Contents

- `data/`: complete point and curve tables, corrections, and statistics
- [`data/published/sources.json`](data/published/sources.json): URLs and SHA-256 hashes of the
  [Bennett–Ghadermarzi tables](https://personal.math.ubc.ca/~bennett/BeGa-data.html)
- `src/`: enumeration, solving, and independent checks
- `patches/`: the PARI/GP patch and Magma cross-check

[`data/points.csv.gz`](data/points.csv.gz) contains signed points with columns `k`, `x`, and `y`.
[`data/curves.jsonl.gz`](data/curves.jsonl.gz) contains a record for every nonzero coefficient in
the range, including empty curves: `points_y_nonnegative` lists pairs up to sign, and
`affine_point_count` counts both signs, counting zero ordinates once. `data/corrected-curves.json`
uses the same nonnegative-ordinate convention. Locally generated comparison snapshots contain rows
of the form `[k, [[x, y], ...]]`, with both signs included.

## Verification

All commands run from the repository root. Python scripts require Python 3.10 or later, use the
standard library, and must run without `-O` so that assertions remain enabled.

To prepare the comparison snapshots, install Poppler’s `pdftotext` and run the first command once.
It downloads the 20 source PDFs and checks their SHA-256 hashes before extracting the rows. If you
already have the audited PDFs, pass `--pdf-dir DIRECTORY` with the filenames from the manifest.

```sh
python3 src/prepare_published.py
python3 verify.py
```

The verifier checks point equations, coefficient coverage, agreement between tables, differences
from the published data, statistics, and the printed-formula counterexamples. This checks the stored
data; it does not repeat the Thue solves.

## Reproduction

Recomputation additionally requires a C++17 compiler supporting `__int128` and PARI/GP 2.19.0 with
the [patch below](#parigp-patch).

### An interval

Prepare the snapshots as described under “[Verification](#verification).” This example covers
$-1000 \le k \le -1$. Bounds are inclusive; $k = 0$ is excluded. For the full range, use disjoint
intervals covering $0 < |k| \le 10^7$.

```sh
mkdir -p work/bin work/example/enumeration
c++ -O3 -std=c++17 src/enumerate_forms.cpp -o work/bin/enumerate
c++ -O3 -std=c++17 src/search_box.cpp -o work/bin/box

work/bin/enumerate -1000 -1 \
  work/example/enumeration/forms.tsv \
  work/example/enumeration/reducible.tsv \
  > work/example/enumeration/enumeration.json

python3 src/check_enumeration.py \
  --enumeration-dir work/example/enumeration \
  --sample-size 1000 --workers 2 \
  --output work/example/enumeration/enumeration-controls.json
python3 src/classify_forms.py \
  --forms work/example/enumeration/forms.tsv \
  --output work/example/enumeration/classification.json

gzip -dc data/published/negative-rows.jsonl.gz > work/negative-rows.jsonl
python3 src/solve_interval.py \
  --enumeration-dir work/example/enumeration \
  --output-dir work/example/solved --workers 2 --gp gp \
  --published-rows work/negative-rows.jsonl

python3 src/verify_interval.py \
  --enumeration-dir work/example/enumeration \
  --result-dir work/example/solved \
  --reducibility-check work/example/enumeration/classification.json \
  --affected-reference data/corrected-curves.json

work/bin/box -1000 -1 100000000 \
  work/example/enumeration/box.tsv \
  > work/example/enumeration/box.json
python3 src/compare_interval.py \
  --result-dir work/example/solved \
  --enumeration-dir work/example/enumeration \
  --published-rows work/negative-rows.jsonl \
  --known-missing data/missing-points.csv
```

Use the positive snapshot for positive intervals. For intervals crossing zero, pass both snapshots
to the solver and concatenate them for the comparison script. Set `--workers` between 1 and 32; add
`--resume` when reusing a previous run’s output directory, with compatible code and inputs.

`--sample-size` controls how many coefficients receive the independent enumeration check; the
example checks every coefficient. The coordinate search is bounded at $x = 10^8$; completeness
relies on enumeration and unconditional Thue solving. Retain the job files until
`verify_interval.py` and `compare_interval.py` pass and the outputs are saved.

### Individual curves

The independent Python solver accepts a list of coefficients:

```sh
python3 src/solve_selected_curves.py -432 -71 \
  --output-dir work/selected --gp gp
```

## PARI/GP patch

The patch replaces `floorr(B0)` with `myfloor(B0)` in `get_B0` (`src/modules/thue.c`), where eleven
equations failed with `precision too low in truncr`. Apply it from the root of the
[official PARI/GP 2.19.0 source tree](https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.0.tar.gz):

```sh
patch -p1 < /path/to/mordell/patches/pari-thue-bound.patch
```

Build PARI according to its installation instructions, run `make test-thue test-kernel`, and supply
the resulting executable with `--gp`.

[The Magma check](patches/magma-check.m), run with Magma V2.29-10 and
`Solutions(Thue(P), 1 : Exact := true)`, agreed with patched PARI on all eleven equations. The
cross-check results are in [`magma-results.txt`](patches/magma-results.txt) and
[`pari-results.txt`](patches/pari-results.txt).

## Algebraic exception

For $k = -432$, the solver uses the algebraic proof in the note to handle a separate PARI/GP
`thue()` failure.
