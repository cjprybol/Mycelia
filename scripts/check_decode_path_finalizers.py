#!/usr/bin/env python3
"""Fail if a decode-path module introduces finalizer, WeakRef, or WeakKeyDict.

The opt5 between-batch GC knob is output-neutral only while the decoder does
not observe collection. A finalizer or weak reference in a decode-path module
would let ``GC.gc()`` between batches change results. None of these tokens
exist on the decode path today (td-3vge item 5). This guard fails when one
shows up in code.

Comments, strings, character literals, and command literals are masked before
the search, so documenting the ban does not trip the guard and a ``#`` inside
a string cannot hide a later call. ``eval`` of a constructed token is out of
scope.

The scanned set is the iterative corrector's per-read decode (the code that
runs between opt5's between-batch ``GC.gc()``) plus the Viterbi modules and
batched kernels that implement it:

- ``src/iterative-assembly.jl``
- ``src/viterbi-next.jl``
- ``src/viterbi-polishing-and-error-correction.jl``
- ``src/rhizomorph/algorithms/path-finding.jl`` (``path_to_sequence``)
- ``src/rhizomorph/algorithms/batched-viterbi-poc.jl``
- ``src/rhizomorph/algorithms/batched-viterbi-kernel.jl``
- ``ext/MyceliaKernelAbstractionsExt.jl``

``src/viterbi*.jl`` and ``src/rhizomorph/algorithms/batched-viterbi*.jl`` are
also globbed, so a new file with those names is scanned without an edit here.
Removing or renaming a required file fails closed until this list is updated.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

TOKENS = re.compile(r"\b(?:finalizer|WeakRef|WeakKeyDict)\b")

# Required files. A missing path is a lint failure so a rename cannot silently
# drop a module out of the guard.
REQUIRED = (
    "src/iterative-assembly.jl",
    "src/viterbi-next.jl",
    "src/viterbi-polishing-and-error-correction.jl",
    "src/rhizomorph/algorithms/path-finding.jl",
    "src/rhizomorph/algorithms/batched-viterbi-poc.jl",
    "src/rhizomorph/algorithms/batched-viterbi-kernel.jl",
    "ext/MyceliaKernelAbstractionsExt.jl",
)

# Extra names a future decoder file is likely to use. Hits are unioned with
# REQUIRED. A glob that matches nothing is fine; a missing REQUIRED path is not.
GLOBS = (
    "src/viterbi*.jl",
    "src/rhizomorph/algorithms/batched-viterbi*.jl",
)


def mask_non_code(text: str) -> str:
    """Replace comments and literals with spaces, preserving newlines."""
    chars = list(text)
    n = len(chars)
    i = 0

    def blank(start: int, end: int) -> None:
        end = min(end, n)
        for j in range(start, end):
            if chars[j] != "\n":
                chars[j] = " "

    while i < n:
        c = chars[i]
        if c == "#":
            start = i
            if i + 1 < n and chars[i + 1] == "=":
                depth = 1
                i += 2
                while i < n and depth:
                    if i + 1 < n and chars[i] == "#" and chars[i + 1] == "=":
                        depth += 1
                        i += 2
                    elif i + 1 < n and chars[i] == "=" and chars[i + 1] == "#":
                        depth -= 1
                        i += 2
                    else:
                        i += 1
            else:
                while i < n and chars[i] != "\n":
                    i += 1
            blank(start, i)
            continue
        if c == '"':
            start = i
            if i + 2 < n and chars[i + 1] == '"' and chars[i + 2] == '"':
                i += 3
                while i < n:
                    if (
                        i + 2 < n
                        and chars[i] == '"'
                        and chars[i + 1] == '"'
                        and chars[i + 2] == '"'
                    ):
                        i += 3
                        break
                    i += 1
            else:
                i += 1
                while i < n:
                    if chars[i] == "\\":
                        i += 2
                        continue
                    if chars[i] == '"' or chars[i] == "\n":
                        if chars[i] == '"':
                            i += 1
                        break
                    i += 1
            blank(start, i)
            continue
        if c == "'":
            # Character literal. Copy it through so '#' does not start a comment.
            i += 1
            if i < n and chars[i] == "\\":
                i += 2
            elif i < n:
                i += 1
            if i < n and chars[i] == "'":
                i += 1
            continue
        if c == "`":
            start = i
            i += 1
            while i < n:
                if chars[i] == "\\":
                    i += 2
                    continue
                if chars[i] == "`":
                    i += 1
                    break
                i += 1
            blank(start, i)
            continue
        i += 1
    return "".join(chars)


def find_tokens(path: Path) -> list[tuple[int, str, str]]:
    """Return (line, token, line_text) hits in code of ``path``."""
    text = path.read_text(encoding="utf-8")
    masked = mask_non_code(text)
    original = text.splitlines()
    hits: list[tuple[int, str, str]] = []
    for lineno, line in enumerate(masked.splitlines(), 1):
        for match in TOKENS.finditer(line):
            shown = original[lineno - 1].strip() if lineno - 1 < len(original) else ""
            hits.append((lineno, match.group(0), shown))
    return hits


def decode_path_files(root: Path) -> tuple[list[Path], list[str]]:
    """Resolve the decode-path set. The second list is missing required paths."""
    missing = [rel for rel in REQUIRED if not (root / rel).is_file()]
    found: dict[str, Path] = {}
    for rel in REQUIRED:
        path = root / rel
        if path.is_file():
            found[rel] = path
    for pattern in GLOBS:
        for path in sorted(root.glob(pattern)):
            if path.is_file():
                rel = path.relative_to(root).as_posix()
                found[rel] = path
    ordered = [found[rel] for rel in sorted(found)]
    return ordered, missing


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Repository root (default: parent of this script). "
        "Ignored when files are passed explicitly.",
    )
    parser.add_argument(
        "files",
        nargs="*",
        type=Path,
        help="Scan only these files (tests). Default: the decode-path set.",
    )
    args = parser.parse_args(argv)

    missing: list[str] = []
    if args.files:
        paths = list(args.files)
    else:
        root = args.root if args.root is not None else Path(__file__).resolve().parents[1]
        paths, missing = decode_path_files(root)

    for path in paths:
        print(path.as_posix())

    errors: list[str] = []
    for rel in missing:
        errors.append(f"missing decode-path file: {rel}")
    for path in paths:
        if not path.is_file():
            errors.append(f"missing decode-path file: {path.as_posix()}")
            continue
        try:
            hits = find_tokens(path)
        except OSError as exc:
            errors.append(f"{path.as_posix()}: {exc}")
            continue
        for lineno, token, shown in hits:
            errors.append(f"{path.as_posix()}:{lineno}: {token}: {shown}")

    if errors:
        print(
            "decode-path GC-neutrality lint: finalizer, WeakRef, and "
            "WeakKeyDict are not allowed in decode-path modules",
            file=sys.stderr,
        )
        for line in errors:
            print(line, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
