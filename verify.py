#!/usr/bin/env python3
"""Check the complete data, published comparison, statistics and printed criteria."""

import subprocess
import sys
from pathlib import Path

from src.verify_data import verify_data


def main():
    root = Path(__file__).resolve().parent
    for sign in ("negative", "positive"):
        if not (root / "data" / "published" / f"{sign}-rows.jsonl.gz").is_file():
            sys.exit(
                "Prepare the comparison snapshots with python3 src/prepare_published.py first."
            )
    verify_data(root)
    subprocess.run(
        [sys.executable, str(root / "src/check_printed_formulas.py")],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    print("Verified exact-arithmetic checks of the printed criteria.")


if __name__ == "__main__":
    main()
