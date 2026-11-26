#!/usr/bin/env python3
"""Verify every ``benchmarking/**/*_test.jl`` file has a ``test/`` entry point.

``Pkg.test()`` executes ``test/runtests.jl``, which only walks ``test/1_*``
through ``test/8_*`` (plus explicit includes under ``test/``). Nothing under
``benchmarking/`` runs in CI unless a wrapper file under ``test/`` includes it.

Heuristic
---------
1. Recursively discover all ``*_test.jl`` files under ``<root>/benchmarking``.
2. Tokenize every ``*.jl`` file under ``<root>/test`` after stripping ``# ...``
   line comments and nested ``#= ... =#`` block comments (while respecting
   Julia string, command, and character literals so comment markers inside
   quotes or ``# include(...)`` examples in header comments are not misparsed).
3. Extract balanced ``include(...)`` call expressions (single-line or
   multi-line, including ``Base.include(...)``), collect their double-quoted
   string literal arguments in order, split them on ``/`` into path segments,
   and mark a ``benchmarking/.../<name>_test.jl`` file as wrapped when an
   ``include(...)`` call's trailing path segments match its repo-relative path
   (for top-level files, mentioning ``"benchmarking"`` followed by the file's
   basename). This handles both ``joinpath(@__DIR__, "..", "..",
   "benchmarking", "<name>")`` and plain ``"../../benchmarking/<name>"``, as
   well as nested subdirectories such as ``benchmarking/sub/x_test.jl``.

Heuristic limits
----------------
- ``include(...)`` arguments constructed via variables, constants, or string
  interpolations containing inner quotes (where ``"benchmarking"`` or path
  segments do not appear directly as string literals inside ``include(...)``)
  are not recognized.
- Only static ``include(...)`` presence under ``test/`` is checked; control
  flow in ``test/runtests.jl`` (such as ``skip_files``, environment-variable
  gates, ``if false``, or un-swept directories like ``test/in_development/``
  and ``test/deprecated/``) is not evaluated.

Merge order with PR #524
------------------------
``ALLOWLIST`` contains ``benchmarking/d5_ont_adaptive_grid_axes_test.jl``,
which is unwrapped on master. PR #524 adds
``test/4_assembly/d5_ont_adaptive_grid_axes_ci_test.jl``, and that file
includes the suite. ``PENDING_TEST_ENTRYPOINTS`` records that pair so this
check accepts either tree:

- wrapper absent (this PR merged, #524 not yet): the suite stays allowlisted
- wrapper present and its ``include(...)`` reaches the suite (#524 merged
  before or after this lint): the suite is wrapped, and that is not a
  stale-allowlist error

A wrapper file that does not include the suite is still stale. Delete the
allowlist entry and the pending mapping together once #524 is on master.
"""

from __future__ import annotations

import argparse
from collections.abc import Collection, Iterator, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
import sys

# Legacy debt: pre-existing unwrapped benchmarking/*_test.jl suites that have
# not yet been wired into Pkg.test(). Do NOT add new files here — every new
# benchmarking/**/*_test.jl file must have a wrapper under test/. Remove entries
# as wrappers are added or files are deleted; this set may only shrink.
ALLOWLIST: frozenset[str] = frozenset(
    {
        "benchmarking/benchmark_kmer_graph_fixture_test.jl",
        "benchmarking/calibration_metrics_test.jl",
        "benchmarking/d5_ont_adaptive_grid_axes_test.jl",
        "benchmarking/indel_benchmark_common_test.jl",
        "benchmarking/rhizomorph_correction_accuracy_benchmark_test.jl",
        "benchmarking/rhizomorph_correction_accuracy_metrics_test.jl",
        "benchmarking/rhizomorph_correction_validation_sweep_test.jl",
        "benchmarking/rhizomorph_correction_validation_sweep_wiring_test.jl",
    }
)

# In-flight test/ wrappers that are not on master yet. Key: benchmarking suite
# that is also in ALLOWLIST. Value: the wrapper path the other PR adds.
# See "Merge order with PR #524" above. Do not add entries for suites that
# already have a wrapper on this branch; those belong off ALLOWLIST.
PENDING_TEST_ENTRYPOINTS: dict[str, str] = {
    "benchmarking/d5_ont_adaptive_grid_axes_test.jl": (
        "test/4_assembly/d5_ont_adaptive_grid_axes_ci_test.jl"
    ),
}


@dataclass(frozen=True)
class AuditResult:
    """Outcome of auditing ``benchmarking/**/*_test.jl`` entry points."""

    wrapped: tuple[str, ...]
    allowlisted_unwrapped: tuple[str, ...]
    new_unwrapped: tuple[str, ...]
    stale_wrapped: tuple[str, ...]
    stale_missing: tuple[str, ...]
    # Allowlisted suites that are wrapped by a recorded in-flight entry point
    # (PENDING_TEST_ENTRYPOINTS). Present so either merge order with that PR
    # passes; they are not errors.
    order_independent_wrapped: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return not self.new_unwrapped and not self.stale_wrapped and not self.stale_missing


def _scan_delimited(source: str, start: int, quote: str, triple: bool) -> tuple[str, int]:
    """Scan a single- or triple-delimited Julia string or command literal."""
    delim = quote * 3 if triple else quote
    step = 3 if triple else 1
    n = len(source)
    i = start + step
    content_start = i
    while i < n:
        if source[i] == "\\" and i + 1 < n:
            i += 2
            continue
        if triple and source[i : i + 3] == delim:
            # Consume extra consecutive quotes so `"""a""""` closes on the last three.
            while i + 3 < n and source[i + 3] == quote:
                i += 1
            return source[content_start:i], i + 3
        if not triple and source[i] == quote:
            return source[content_start:i], i + 1
        i += 1
    return source[content_start:n], n


def tokenize_julia(source: str) -> Iterator[tuple[str, str]]:
    """Yield simplified tokens from Julia source with all comments stripped."""
    n = len(source)
    i = 0
    prev_can_adjoint = False

    while i < n:
        ch = source[i]

        # Nested block comment: #= ... =#
        if ch == "#" and source[i : i + 2] == "#=":
            depth = 1
            i += 2
            while i < n and depth > 0:
                if source[i : i + 2] == "#=":
                    depth += 1
                    i += 2
                elif source[i : i + 2] == "=#":
                    depth -= 1
                    i += 2
                else:
                    i += 1
            prev_can_adjoint = False
            continue

        # Single-line comment: # ...
        if ch == "#":
            newline = source.find("\n", i + 1)
            i = n if newline == -1 else newline + 1
            prev_can_adjoint = False
            continue

        # Whitespace
        if ch.isspace():
            i += 1
            prev_can_adjoint = False
            continue

        # String literal: """...""" or "..."
        if ch == '"':
            triple = source[i : i + 3] == '"""'
            text, i = _scan_delimited(source, i, '"', triple=triple)
            yield ("STRING", text)
            prev_can_adjoint = False
            continue

        # Command literal: ```...``` or `...`
        if ch == "`":
            triple = source[i : i + 3] == "```"
            _, i = _scan_delimited(source, i, "`", triple=triple)
            prev_can_adjoint = False
            continue

        # Adjoint operator `'` vs character literal `'c'`
        if ch == "'":
            if prev_can_adjoint:
                i += 1
                prev_can_adjoint = True
            else:
                i += 1
                while i < n and source[i] not in ("'", "\n"):
                    if source[i] == "\\" and i + 1 < n:
                        i += 2
                    else:
                        i += 1
                if i < n and source[i] == "'":
                    i += 1
                prev_can_adjoint = False
            continue

        if ch == "(":
            yield ("LPAREN", ch)
            i += 1
            prev_can_adjoint = False
            continue

        if ch == ")":
            yield ("RPAREN", ch)
            i += 1
            prev_can_adjoint = True
            continue

        if ch in ("]", "}"):
            yield ("OTHER", ch)
            i += 1
            prev_can_adjoint = True
            continue

        if ch.isalpha() or ch == "_":
            start = i
            i += 1
            while i < n and (source[i].isalnum() or source[i] in ("_", "!")):
                i += 1
            yield ("IDENT", source[start:i])
            prev_can_adjoint = True
            continue

        if ch.isdigit():
            i += 1
            while i < n and (source[i].isalnum() or source[i] == "_"):
                i += 1
            prev_can_adjoint = True
            continue

        yield ("OTHER", ch)
        i += 1
        prev_can_adjoint = False


def extract_include_segment_tuples(source: str) -> list[tuple[str, ...]]:
    """Extract normalized path-segment tuples from ``include(...)`` calls."""
    tokens = list(tokenize_julia(source))
    results: list[tuple[str, ...]] = []
    n = len(tokens)
    i = 0

    while i < n:
        kind, value = tokens[i]
        if kind == "IDENT" and value == "include" and i + 1 < n and tokens[i + 1][0] == "LPAREN":
            depth = 1
            j = i + 2
            string_literals: list[str] = []
            while j < n and depth > 0:
                tkind, tval = tokens[j]
                if tkind == "LPAREN":
                    depth += 1
                elif tkind == "RPAREN":
                    depth -= 1
                elif tkind == "STRING":
                    string_literals.append(tval)
                j += 1

            segments: list[str] = []
            for literal in string_literals:
                for part in literal.replace("\\", "/").split("/"):
                    if part and part != ".":
                        segments.append(part)
            if segments:
                results.append(tuple(segments))
            i = j
            continue
        i += 1

    return results


def discover_benchmark_tests(root: Path) -> set[str]:
    """Return repo-relative POSIX paths of all ``benchmarking/**/*_test.jl`` files."""
    bench_dir = root / "benchmarking"
    if not bench_dir.is_dir():
        return set()
    discovered: set[str] = set()
    for path in bench_dir.rglob("*_test.jl"):
        if path.is_file() and ".ipynb_checkpoints" not in path.parts:
            discovered.add(path.relative_to(root).as_posix())
    return discovered


def collect_test_include_segments(root: Path) -> list[tuple[str, ...]]:
    """Collect ``include(...)`` path-segment tuples across ``test/**/*.jl``."""
    test_dir = root / "test"
    if not test_dir.is_dir():
        return []
    all_segments: list[tuple[str, ...]] = []
    for jl_file in sorted(test_dir.rglob("*.jl")):
        if not jl_file.is_file() or ".ipynb_checkpoints" in jl_file.parts:
            continue
        source = jl_file.read_text(encoding="utf-8", errors="replace")
        all_segments.extend(extract_include_segment_tuples(source))
    return all_segments


def is_wrapped(rel_posix_path: str, include_segments: Sequence[tuple[str, ...]]) -> bool:
    """Return True if any ``test/`` ``include(...)`` call ends with ``rel_posix_path``."""
    target = PurePosixPath(rel_posix_path).parts
    target_len = len(target)
    return any(
        len(segs) >= target_len and segs[-target_len:] == target
        for segs in include_segments
    )


def pending_wrapper_covers(
    root: Path,
    rel_posix_path: str,
    pending_entrypoints: Mapping[str, str] | None = None,
) -> bool:
    """Return True when a recorded in-flight wrapper includes ``rel_posix_path``.

    The wrapper path has to be the one named in ``PENDING_TEST_ENTRYPOINTS``
    (or ``pending_entrypoints``) and that file's ``include(...)`` calls have to
    reach the suite. A same-named file that does not include it does not count.
    """
    mapping = (
        PENDING_TEST_ENTRYPOINTS
        if pending_entrypoints is None
        else dict(pending_entrypoints)
    )
    wrapper_rel = mapping.get(rel_posix_path)
    if wrapper_rel is None:
        return False
    wrapper_path = root / wrapper_rel
    if not wrapper_path.is_file():
        return False
    source = wrapper_path.read_text(encoding="utf-8", errors="replace")
    return is_wrapped(rel_posix_path, extract_include_segment_tuples(source))


def audit_entrypoints(
    root: Path,
    allowlist: Collection[str] | None = None,
    pending_entrypoints: Mapping[str, str] | None = None,
) -> AuditResult:
    """Audit ``benchmarking/**/*_test.jl`` files under ``root`` against ``test/``."""
    active_allowlist = frozenset(ALLOWLIST if allowlist is None else allowlist)
    discovered = discover_benchmark_tests(root)
    include_segments = collect_test_include_segments(root)

    wrapped = sorted(p for p in discovered if is_wrapped(p, include_segments))
    wrapped_set = set(wrapped)
    unwrapped_set = discovered - wrapped_set

    # Allowlisted suites whose recorded #524-style wrapper is actually present
    # are wrapped on purpose. They are not stale, so either merge order passes.
    order_independent = {
        path
        for path in active_allowlist & wrapped_set
        if pending_wrapper_covers(root, path, pending_entrypoints)
    }
    allowlisted_unwrapped = sorted(unwrapped_set & active_allowlist)
    new_unwrapped = sorted(unwrapped_set - active_allowlist)
    stale_wrapped = sorted((active_allowlist & wrapped_set) - order_independent)
    stale_missing = sorted(active_allowlist - discovered)

    return AuditResult(
        wrapped=tuple(wrapped),
        allowlisted_unwrapped=tuple(allowlisted_unwrapped),
        new_unwrapped=tuple(new_unwrapped),
        stale_wrapped=tuple(stale_wrapped),
        stale_missing=tuple(stale_missing),
        order_independent_wrapped=tuple(sorted(order_independent)),
    )


def _format_wrapper_example(rel_posix_path: str) -> str:
    parts = PurePosixPath(rel_posix_path).parts
    quoted_parts = ", ".join(f'"{part}"' for part in parts)
    return (
        "  module MyBenchmarkSuiteTests\n"
        f'  include(joinpath(@__DIR__, "..", "..", {quoted_parts}))\n'
        "  end"
    )


def main(
    argv: Sequence[str] | None = None,
    *,
    allowlist: Collection[str] | None = None,
) -> int:
    parser = argparse.ArgumentParser(
        description="Check that benchmarking/**/*_test.jl files have a test/ entry point."
    )
    default_root = Path(__file__).resolve().parent.parent
    parser.add_argument(
        "--root",
        type=Path,
        default=default_root,
        help="Repository root directory to inspect (default: %(default)s).",
    )
    args = parser.parse_args(argv)
    root = args.root.resolve()

    if not root.is_dir():
        print(f"ERROR: Repository root does not exist: {root}", file=sys.stderr)
        return 1

    result = audit_entrypoints(root, allowlist=allowlist)

    if result.new_unwrapped:
        print(
            "ERROR: Benchmarking test file(s) have no include(...) entry point under test/:",
            file=sys.stderr,
        )
        for path in result.new_unwrapped:
            print(f"  - {path}", file=sys.stderr)
        print(
            "\nPkg.test() does not walk benchmarking/. Add a thin module wrapper under\n"
            "test/ (for example in test/4_assembly/<name>_test.jl, following\n"
            "test/4_assembly/ont_k_sweep_shards_driver_test.jl):\n",
            file=sys.stderr,
        )
        print(_format_wrapper_example(result.new_unwrapped[0]), file=sys.stderr)

    if result.stale_wrapped or result.stale_missing:
        if result.new_unwrapped:
            print(file=sys.stderr)
        print(
            "ERROR: Stale entry/entries in ALLOWLIST in ci/check_benchmark_test_entrypoints.py "
            "(the allowlist may only shrink):",
            file=sys.stderr,
        )
        for path in result.stale_wrapped:
            print(f"  - {path} (now wrapped under test/; remove from ALLOWLIST)", file=sys.stderr)
        for path in result.stale_missing:
            print(f"  - {path} (file no longer exists; remove from ALLOWLIST)", file=sys.stderr)

    if not result.ok:
        return 1

    if result.order_independent_wrapped:
        print(
            "NOTE: allowlisted suite(s) are wrapped by a recorded in-flight test/ "
            "entry point (PENDING_TEST_ENTRYPOINTS). This is accepted so the "
            "allowlist passes whether or not that wrapper is in the tree. "
            "Remove the suite from ALLOWLIST and PENDING_TEST_ENTRYPOINTS once "
            "that entry point is on master:"
        )
        for path in result.order_independent_wrapped:
            wrapper = PENDING_TEST_ENTRYPOINTS[path]
            print(f"  - {path} via {wrapper}")

    total = len(result.wrapped) + len(result.allowlisted_unwrapped)
    print(
        f"OK: {total} benchmarking test file(s) checked "
        f"({len(result.wrapped)} wrapped under test/, "
        f"{len(result.allowlisted_unwrapped)} allowlisted legacy debt)."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
