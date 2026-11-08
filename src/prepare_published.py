"""Rebuild the comparison snapshots from the audited Bennett–Ghadermarzi PDFs."""

import argparse
import gzip
import json
import re
import shutil
import subprocess
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

from common import sha256_file

PUBLISHED = Path(__file__).resolve().parents[1] / "data" / "published"
HEADER = re.compile(r"^[ \t]*(-?\d+)[ \t]+(\d+)\s+(?=\[)", re.MULTILINE)
POINT = re.compile(r"\[\s*(-?\d+)\s*,\s*(-?\d+)\s*\]")


def write_rows(pdf, table, output):
    text = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"]).decode("utf-8")
    prefix, *fields = HEADER.split(text.replace("\f", "\n"))
    if prefix.strip() or not fields:
        raise ValueError(f"Unexpected table layout in {pdf.name}.")
    previous = table["lower_exclusive"]
    for index in range(0, len(fields), 3):
        k_text, count, body = fields[index : index + 3]
        k = int(k_text)
        magnitude = k * table["sign"]
        points = [(int(x), int(y)) for x, y in POINT.findall(body)]
        if not previous < magnitude <= table["upper_inclusive"]:
            raise ValueError(f"Unexpected coefficient {k} in {pdf.name}.")
        if POINT.sub("", body).strip() or len(points) != int(count):
            raise ValueError(f"Invalid point list at k = {k} in {pdf.name}.")
        output.write(json.dumps([k, points], separators=(",", ":")) + "\n")
        previous = magnitude


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf-dir", type=Path, help="Read the named PDFs from this directory.")
    args = parser.parse_args()
    if shutil.which("pdftotext") is None:
        parser.exit(1, "Install Poppler’s pdftotext before preparing the snapshots.\n")
    tables = json.loads((PUBLISHED / "sources.json").read_text())["tables"]
    with TemporaryDirectory(dir=PUBLISHED) as staging:
        stage = Path(staging)
        for sign, name in [(-1, "negative"), (1, "positive")]:
            with gzip.open(
                stage / f"{name}-rows.jsonl.gz",
                "wt",
                encoding="utf-8",
                newline="\n",
                compresslevel=6,
            ) as output:
                for table in tables:
                    if table["sign"] != sign:
                        continue
                    pdf = (args.pdf_dir or stage) / table["file"]
                    print(f"Reading {pdf.name}…", flush=True)
                    if args.pdf_dir is None:
                        with urlopen(table["url"], timeout=60) as response, pdf.open("wb") as dest:
                            shutil.copyfileobj(response, dest)
                    if sha256_file(pdf) != table["sha256"]:
                        raise ValueError(f"{pdf.name} differs from the audited PDF.")
                    write_rows(pdf, table, output)
        for snapshot in stage.glob("*-rows.jsonl.gz"):
            snapshot.replace(PUBLISHED / snapshot.name)


if __name__ == "__main__":
    main()
