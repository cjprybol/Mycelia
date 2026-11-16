#!/usr/bin/env python3
"""Unit tests for ``ci/check_benchmark_test_entrypoints.py``."""

from __future__ import annotations

import contextlib
import io
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import check_benchmark_test_entrypoints as lint  # noqa: E402


def _write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


class CheckBenchmarkTestEntrypointsTests(unittest.TestCase):
    def _run_main(
        self,
        root: Path,
        allowlist: frozenset[str] = frozenset(),
    ) -> tuple[int, str, str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with (
            mock.patch.object(lint, "ALLOWLIST", allowlist),
            contextlib.redirect_stdout(stdout),
            contextlib.redirect_stderr(stderr),
        ):
            rc = lint.main(["--root", str(root)])
        return rc, stdout.getvalue(), stderr.getvalue()

    def test_wrapped_file_joinpath_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root / "benchmarking" / "alpha_test.jl", "# suite\n")
            _write(
                root / "test" / "4_assembly" / "alpha_driver_test.jl",
                'module AlphaDriverTests\n'
                'include(joinpath(@__DIR__, "..", "..", "benchmarking", "alpha_test.jl"))\n'
                'end\n',
            )

            rc, out, err = self._run_main(root)
            self.assertEqual(rc, 0, msg=err)
            self.assertEqual(err, "")
            self.assertIn("1 wrapped under test/", out)

    def test_wrapped_file_plain_relative_and_multiline_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root / "benchmarking" / "plain_test.jl", "# suite\n")
            _write(root / "benchmarking" / "multiline_test.jl", "# suite\n")
            _write(
                root / "test" / "4_assembly" / "wrappers_test.jl",
                'module PlainWrapper\n'
                'include("../../benchmarking/plain_test.jl")\n'
                'end\n'
                'module MultilineWrapper\n'
                'include(\n'
                '    joinpath(\n'
                '        @__DIR__,\n'
                '        "..",\n'
                '        "..",\n'
                '        "benchmarking",\n'
                '        "multiline_test.jl",\n'
                '    ),\n'
                ')\n'
                'end\n',
            )

            rc, _, err = self._run_main(root)
            self.assertEqual(rc, 0, msg=err)

    def test_new_unwrapped_file_fails_and_is_named(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root / "benchmarking" / "unwrapped_new_suite_test.jl", "# suite\n")
            (root / "test").mkdir()

            rc, _, err = self._run_main(root)
            self.assertEqual(rc, 1)
            self.assertIn("benchmarking/unwrapped_new_suite_test.jl", err)
            self.assertIn("module ", err)
            self.assertIn('include(joinpath(@__DIR__, "..", "..", "benchmarking", "unwrapped_new_suite_test.jl"))', err)

    def test_allowlisted_unwrapped_file_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root / "benchmarking" / "legacy_unwrapped_test.jl", "# legacy\n")
            (root / "test").mkdir()

            rc, out, err = self._run_main(
                root,
                allowlist=frozenset({"benchmarking/legacy_unwrapped_test.jl"}),
            )
            self.assertEqual(rc, 0, msg=err)
            self.assertEqual(err, "")
            self.assertIn("1 allowlisted legacy debt", out)

    def test_stale_allowlist_entry_now_wrapped_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root / "benchmarking" / "now_wrapped_test.jl", "# suite\n")
            _write(
                root / "test" / "4_assembly" / "now_wrapped_driver_test.jl",
                'module NowWrappedTests\n'
                'include(joinpath(@__DIR__, "..", "..", "benchmarking", "now_wrapped_test.jl"))\n'
                'end\n',
            )

            rc, _, err = self._run_main(
                root,
                allowlist=frozenset({"benchmarking/now_wrapped_test.jl"}),
            )
            self.assertEqual(rc, 1)
            self.assertIn("benchmarking/now_wrapped_test.jl", err)
            self.assertIn("ALLOWLIST", err)

    def test_stale_allowlist_entry_file_deleted_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "benchmarking").mkdir()
            (root / "test").mkdir()

            rc, _, err = self._run_main(
                root,
                allowlist=frozenset({"benchmarking/deleted_suite_test.jl"}),
            )
            self.assertEqual(rc, 1)
            self.assertIn("benchmarking/deleted_suite_test.jl", err)
            self.assertIn("ALLOWLIST", err)

    def test_nested_benchmark_test_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root / "benchmarking" / "sub" / "x_test.jl", "# nested suite\n")
            (root / "test").mkdir()

            # Unwrapped nested _test.jl fails and is named with its subdirectory path.
            rc, _, err = self._run_main(root)
            self.assertEqual(rc, 1)
            self.assertIn("benchmarking/sub/x_test.jl", err)

            # Wrapping a sibling benchmarking/other/x_test.jl does NOT satisfy sub/x_test.jl.
            _write(root / "benchmarking" / "other" / "x_test.jl", "# sibling\n")
            _write(
                root / "test" / "4_assembly" / "other_x_driver_test.jl",
                'module OtherXTests\n'
                'include(joinpath(@__DIR__, "..", "..", "benchmarking", "other", "x_test.jl"))\n'
                'end\n',
            )
            rc2, _, err2 = self._run_main(root)
            self.assertEqual(rc2, 1)
            self.assertIn("benchmarking/sub/x_test.jl", err2)
            self.assertNotIn("benchmarking/other/x_test.jl", err2)

            # Wrapping benchmarking/sub/x_test.jl makes both pass.
            _write(
                root / "test" / "4_assembly" / "sub_x_driver_test.jl",
                'module SubXTests\n'
                'include(joinpath(@__DIR__, "..", "..", "benchmarking", "sub", "x_test.jl"))\n'
                'end\n',
            )
            rc3, _, err3 = self._run_main(root)
            self.assertEqual(rc3, 0, msg=err3)

    def test_commented_out_include_does_not_count_as_entrypoint(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _write(root / "benchmarking" / "commented_test.jl", "# suite\n")
            _write(
                root / "test" / "4_assembly" / "commented_driver_test.jl",
                '# Run: julia -e \'include("../../benchmarking/commented_test.jl")\'\n'
                '# include(joinpath(@__DIR__, "..", "..", "benchmarking", "commented_test.jl"))\n'
                '    # include("../../benchmarking/commented_test.jl")\n'
                '#=\n'
                'include(joinpath(@__DIR__, "..", "..", "benchmarking", "commented_test.jl"))\n'
                '=#\n'
                'const EXAMPLE = """\n'
                'include(joinpath(@__DIR__, "..", "..", "benchmarking", "commented_test.jl"))\n'
                '"""\n',
            )

            rc, _, err = self._run_main(root)
            self.assertEqual(rc, 1)
            self.assertIn("benchmarking/commented_test.jl", err)

    def test_cli_subprocess_with_root_argument(self) -> None:
        script = Path(__file__).resolve().parent / "check_benchmark_test_entrypoints.py"
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for allowlisted in lint.ALLOWLIST:
                _write(root / allowlisted, "# legacy\n")
            _write(root / "benchmarking" / "brand_new_test.jl", "# new\n")
            (root / "test").mkdir()

            proc_fail = subprocess.run(
                [sys.executable, "-B", str(script), "--root", str(root)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_fail.returncode, 1)
            self.assertIn("benchmarking/brand_new_test.jl", proc_fail.stderr)

            _write(
                root / "test" / "brand_new_wrapper_test.jl",
                'include("../benchmarking/brand_new_test.jl")\n',
            )
            proc_ok = subprocess.run(
                [sys.executable, "-B", str(script), "--root", str(root)],
                capture_output=True,
                text=True,
                check=False,
            )
            self.assertEqual(proc_ok.returncode, 0, msg=proc_ok.stderr)


if __name__ == "__main__":
    unittest.main()
