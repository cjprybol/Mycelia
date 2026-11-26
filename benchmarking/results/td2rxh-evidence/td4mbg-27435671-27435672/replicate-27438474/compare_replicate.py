#!/usr/bin/env python3
"""Amendment-6 diagnostic: BASE replicate vs original BASE, same G0 field rule.

Usage: compare_replicate.py <original BASE shard dir> <replicate shard dir> <seed>...
Both directories hold seed<N>.csv files (shards/base/ and replicate-27438474/shards/).
Compares every column except label/wall_seconds/wall_seconds_per_read as exact
strings, per (cell, seed, arm). Prints a summary and every differing field.

Refuses a short argument list, a missing directory or shard, a non-canonical
seed argument, a row whose seed field is not that argument, a run label that is
not td2rxh@<12 hex>+manifest<12 hex>, and two sides that do not share one label.
A BASE replicate that names a different commit is not a same-commit comparison.
"""
import csv
import re
import sys
from pathlib import Path

EXCLUDED = {"label", "wall_seconds", "wall_seconds_per_read"}
LABEL_RE = re.compile(r"^td2rxh@[0-9a-f]{12}\+manifest[0-9a-f]{12}$")
REQUIRED = ("cell", "seed", "arm", "label")


def rows(path: Path) -> dict:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames or []
        missing = [c for c in REQUIRED if c not in fields]
        if missing:
            raise SystemExit(f"{path}: missing columns {missing}")
        loaded = {}
        for record in reader:
            key = (record["cell"], record["seed"], record["arm"])
            if key in loaded:
                raise SystemExit(f"{path}: duplicate row {key}")
            loaded[key] = record
        return loaded


def _one_label(path: Path, table: dict, seed: str) -> str:
    if not table:
        raise SystemExit(f"{path}: no data rows")
    labels = set()
    for key, record in table.items():
        if key[1] != seed:
            raise SystemExit(
                f"{path}: row seed {key[1]!r} does not match requested seed {seed}"
            )
        label = record["label"]
        if LABEL_RE.fullmatch(label) is None:
            raise SystemExit(
                f"{path}: run label {label!r} is not td2rxh@<12 hex>+manifest<12 hex>"
            )
        labels.add(label)
    if len(labels) != 1:
        raise SystemExit(f"{path}: rows carry {len(labels)} run labels")
    return next(iter(labels))


def main() -> None:
    if len(sys.argv) < 4:
        raise SystemExit(
            "usage: compare_replicate.py <original BASE shard dir> "
            "<replicate shard dir> <seed>..."
        )
    orig_dir, rep_dir = Path(sys.argv[1]), Path(sys.argv[2])
    for directory in (orig_dir, rep_dir):
        if not directory.is_dir():
            raise SystemExit(f"not a directory: {directory}")
    seeds = sys.argv[3:]
    for seed in seeds:
        if not seed.isdigit() or str(int(seed)) != seed:
            raise SystemExit(f"seed {seed!r} is not a canonical integer")
    total_rows = differing_rows = 0
    labels = set()
    for seed in seeds:
        orig_path = orig_dir / f"seed{seed}.csv"
        rep_path = rep_dir / f"seed{seed}.csv"
        for path in (orig_path, rep_path):
            if not path.is_file():
                raise SystemExit(f"missing shard: {path}")
        orig = rows(orig_path)
        rep = rows(rep_path)
        orig_label = _one_label(orig_path, orig, seed)
        rep_label = _one_label(rep_path, rep, seed)
        if orig_label != rep_label:
            raise SystemExit(
                f"seed {seed}: labels differ ({orig_label} vs {rep_label})"
            )
        labels.add(orig_label)
        if set(orig) != set(rep):
            raise SystemExit(f"seed {seed}: row keys differ")
        for key in sorted(orig):
            total_rows += 1
            cols = [c for c in orig[key] if c not in EXCLUDED]
            diffs = [(c, orig[key][c], rep[key][c]) for c in cols if orig[key][c] != rep[key][c]]
            if diffs:
                differing_rows += 1
            status = "DIFF" if diffs else "same"
            print(f"{status} {key[0]} seed={key[1]} {key[2]}" +
                  "".join(f" | {c}: {a} -> {b}" for c, a, b in diffs))
    print(f"labels: {sorted(labels)}")
    print(f"rows compared={total_rows} differing={differing_rows}")


if __name__ == "__main__":
    main()
