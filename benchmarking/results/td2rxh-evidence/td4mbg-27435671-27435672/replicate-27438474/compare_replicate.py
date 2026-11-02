#!/usr/bin/env python3
"""Amendment-6 diagnostic: BASE replicate vs original BASE, same G0 field rule.

Usage: compare_rep.py <original BASE shard dir> <replicate shard dir> <seed>...
Compares every column except label/wall_seconds/wall_seconds_per_read as exact
strings, per (cell, seed, arm). Prints a summary and every differing field.
"""
import csv
import sys
from pathlib import Path

EXCLUDED = {"label", "wall_seconds", "wall_seconds_per_read"}


def rows(path: Path) -> dict:
    with open(path, newline="") as fh:
        return {(r["cell"], r["seed"], r["arm"]): r for r in csv.DictReader(fh)}


def main() -> None:
    orig_dir, rep_dir = Path(sys.argv[1]), Path(sys.argv[2])
    total_rows = differing_rows = 0
    labels = set()
    for seed in sys.argv[3:]:
        orig = rows(orig_dir / f"seed{seed}.csv")
        rep = rows(rep_dir / f"rep_seed{seed}.csv")
        labels |= {r["label"] for r in orig.values()} | {r["label"] for r in rep.values()}
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
