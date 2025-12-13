# HPC HPC-Backed CI Plan (Phased)

Phased plan for integrating HPC-driven CI while keeping GitHub Actions as-is. Mirrors the supplied proposal and adds repo-specific context for Mycelia.

## 0. Goals & Invariants
- Preserve existing GitHub Actions CI: lightweight unit tests, Codecov upload, existing README badge remains the baseline signal.
- Add HPC-only coverage, extended correctness tests, and benchmarks; surface results back to GitHub with extra Codecov flag(s), badges, and optional PR/status checks.
- Assume HPC harness can emit coverage (`lcov.info`), benchmark JSON/HTML, and correctness summaries (pass/fail counts, regression info).

## 1. Common Building Blocks (All Phases)
- **Single HPC driver**: `ci/hpc/run_hpc_ci.sh` (or `.jl` wrapper) that:
  - Checks out the intended commit using `HPC_CI_COMMIT` (SHA) and `HPC_CI_BRANCH` (optional for Codecov) via `git fetch origin "$HPC_CI_BRANCH"` then `git checkout "$HPC_CI_COMMIT"`.
  - Runs extended tests with coverage, e.g. `julia --project=. --code-coverage=user ci/run_extended_tests.jl --mode=all --hpc`.
  - Runs benchmarks, e.g. `julia --project=. benchmarking/benchmark_runner.jl medium --hpc` or `benchmarking/run_all_benchmarks.sh`.
  - Writes standardized artifacts under `hpc-ci/artifacts/`:
    - `coverage/lcov.info`
    - `tests/summary.json` (overall pass/fail counts, failing suites)
    - `benchmarks/summary.json` (key metrics per benchmark, averages, regression flags)
    - `meta.json` (commit SHA, branch, Julia version, HPC machine name, date, SLURM job ID, etc.)
- Reuse this script for: manual/cron HPC runs (Phase 1), HPC → GitHub data push (Phase 2), and Jacamar/GitLab runner jobs (Phase 3).
- Repo fit: keep Julia imports centralized (per `AGENTS.md`), prefer `joinpath`, and place new scripts under `ci/hpc/` so they auto-load only when invoked.

## 2. Phase 1 – Immediate Easy Win (Manual/Cron HPC → Codecov Flag)
- **Goal**: Upload HPC coverage to Codecov with a distinct flag/badge; run extended tests/benchmarks out-of-band.
- **Changes**: No GitHub workflow edits beyond adding a badge later. Extend `run_hpc_ci.sh` to upload coverage; optionally keep JSON artifacts on-cluster.
- **Steps**:
  1. Add a Codecov token to the HPC environment (`CODECOV_TOKEN`) stored securely.
  2. After coverage generation, upload with a flag:
     ```bash
     codecov -t "$CODECOV_TOKEN" -F hpc-extended -f hpc-ci/artifacts/coverage/lcov.info -C "$HPC_CI_COMMIT" -B "$HPC_CI_BRANCH"
     ```
     - `-F hpc-extended` defines a dedicated Codecov flag distinct from GitHub CI.
  3. Add an HPC coverage badge next to the existing one, e.g.:
     ```
     [![HPC Coverage](https://codecov.io/gh/ORG/REPO/branch/main/graph/badge.svg?flag=hpc-extended)](...)
     ```
  4. Trigger via manual runs or cron on the HPC login node (e.g., nightly on `main`, release tags).
- **Outcome**: Immediate visibility of HPC-only coverage in Codecov without touching GitHub workflows.

## 3. Phase 2 – Robust Integration Without GitHub Self-Hosted Runner
- **Goal**: Make HPC results first-class on GitHub (badges/docs/status) without registering an HPC node as a GitHub runner.
- **HPC side (push summaries)**:
  1. Create a dedicated branch:
     ```bash
     git checkout --orphan hpc-results
     rm -rf .
     echo "# HPC CI results branch" > README.md
     git add README.md
     git commit -m "Initialize hpc-results branch"
     git push origin hpc-results
     ```
  2. Create a deploy key or PAT scoped to this repo (ideally only to `hpc-results`).
  3. In `run_hpc_ci.sh`, after artifact creation, push small JSON summaries:
     ```bash
     WORKDIR=$(mktemp -d)
     cd "$WORKDIR"
     git clone --single-branch --branch hpc-results git@github.com:ORG/REPO.git .
     mkdir -p "$HPC_CI_COMMIT"
     cp /path/to/hpc-ci/artifacts/tests/summary.json "$HPC_CI_COMMIT/tests.json"
     cp /path/to/hpc-ci/artifacts/benchmarks/summary.json "$HPC_CI_COMMIT/benchmarks.json"
     cp /path/to/hpc-ci/artifacts/meta.json "$HPC_CI_COMMIT/meta.json"
     jq '.' "$HPC_CI_COMMIT/tests.json" > latest-tests.json
     jq '.' "$HPC_CI_COMMIT/benchmarks.json" > latest-benchmarks.json
     git add "$HPC_CI_COMMIT" latest-*.json
     git commit -m "HPC results for $HPC_CI_COMMIT"
     git push origin hpc-results
     ```
     - Only push small structured JSON; keep large HTML/CSV logs on-cluster.
- **GitHub side (ingest/publish)**:
  - Add `.github/workflows/hpc-results.yml` triggered on `push` to `hpc-results`; copy `latest-*.json` to `docs/hpc/` or serve via GitHub Pages from `hpc-results`.
  - Use Shields dynamic badges pointed at these JSON endpoints, e.g.:
    ```
    ![HPC Tests](https://img.shields.io/endpoint?url=https://<gh-pages>/latest-tests.json&label=HPC%20tests)
    ![HPC Benchmarks](https://img.shields.io/endpoint?url=https://<gh-pages>/latest-benchmarks.json&label=HPC%20bench)
    ```
  - Coverage still via the Codecov `hpc-extended` flag badge from Phase 1.
- **Optional GitHub checks/comments**:
  - HPC script can post commit statuses or PR comments (`context: hpc-tests`, `hpc-benchmarks`) pointing to internal logs.
  - Alternatively, a GitHub Action on `hpc-results` can map `meta.json` to the relevant commit/PR and comment with findings.
- **Outcome**: Automated HPC data flow into GitHub badges/docs without any self-hosted runner.

## 4. Phase 3 – Long-Term: Jacamar CI + GitLab Runners on HPC
- **Goal**: Run HPC CI via facility-provided Jacamar-backed GitLab runners while surfacing results to GitHub.
- **Topology**:
  - GitHub stays canonical for issues/PRs/badges; maintain a GitLab mirror (or GitLab primary with GitHub mirror) per facility support.
  - Jacamar-backed GitLab runners submit jobs to the scheduler; reuse `ci/hpc/run_hpc_ci.sh` on compute nodes.
- **`.gitlab-ci.yml` sketch**:
  ```yaml
  stages: [hpc-test, hpc-bench]
  variables: { GIT_SUBMODULE_STRATEGY: recursive }
  hpc-tests:
    stage: hpc-test
    tags: [jacamar]
    script:
      - export HPC_CI_COMMIT="$CI_COMMIT_SHA"
      - export HPC_CI_BRANCH="$CI_COMMIT_REF_NAME"
      - ci/hpc/run_hpc_ci.sh --tests-only
    artifacts:
      when: always
      paths:
        - hpc-ci/artifacts/tests
        - hpc-ci/artifacts/coverage
        - hpc-ci/artifacts/meta.json
  hpc-benchmarks:
    stage: hpc-bench
    needs: [hpc-tests]
    tags: [jacamar]
    script:
      - export HPC_CI_COMMIT="$CI_COMMIT_SHA"
      - ci/hpc/run_hpc_ci.sh --benchmarks-only
    artifacts:
      when: always
      paths:
        - hpc-ci/artifacts/benchmarks
        - hpc-ci/artifacts/meta.json
  publish-hpc-results:
    stage: hpc-bench
    needs: [hpc-tests, hpc-benchmarks]
    tags: [jacamar]
    script:
      - codecov -t "$CODECOV_TOKEN" -F jacamar-hpc -f hpc-ci/artifacts/coverage/lcov.info -C "$CI_COMMIT_SHA" -B "$CI_COMMIT_REF_NAME"
      - ./ci/hpc/push_to_hpc_results_branch.sh
  ```
  - `tags: [jacamar]` targets Jacamar runners. `push_to_hpc_results_branch.sh` mirrors Phase 2 logic, now running inside GitLab CI.
- **Sync back to GitHub**:
  - Codecov as the coverage hub with flags `ci-lite` (GH Actions), `hpc-extended` (Phase 1/2), `jacamar-hpc` (Phase 3). Add badges as desired.
  - Reuse `hpc-results` branch + GitHub Actions for badges/docs; GitLab job pushes summaries using a GitHub PAT stored as a GitLab variable.
  - Optional GitLab job hits GitHub Statuses API (`context: jacamar/hpc-tests`) with links to GitLab pipeline pages.
- **Outcome**: Fully automated HPC CI via facility-supported runners, results visible on GitHub.

## 5. Data Model & Provenance (All Phases)
- **Coverage**: Current CI coverage remains default; add flags `hpc-extended` and `jacamar-hpc`.
- **Test summaries** (`tests/summary.json` example):
  ```json
  {"commit":"<sha>","status":"passed","suites":{"unit":{"passed":123,"failed":0},"integration":{"passed":45,"failed":0},"tutorials":{"passed":10,"failed":1}}}
  ```
- **Benchmark summaries** (`benchmarks/summary.json` example):
  ```json
  {"commit":"<sha>","benchmarks":{"solver/3d":{"median_ms":123.4,"delta_vs_baseline_pct":-2.1},"solver/2d":{"median_ms":56.7,"delta_vs_baseline_pct":0.5}}}
  ```
  - A single “worst regression” number can feed a README badge.
- **Provenance (`meta.json`)**: commit SHA/branch/tag, Julia version, compiler flags/feature toggles, HPC system/partition, walltime/node count, scheduler job ID, timestamp.

## 6. Security & Policy Notes
- **Tokens**: Codecov token stays on HPC/Jacamar runners; GitHub deploy key/PAT scoped to a single repo/branch; GitLab variables masked/protected.
- **Size control**: Push only small JSON to GitHub; keep large logs/artifacts in facility storage or object store.
- **Staleness**: Include timestamps + commit SHA in `latest` JSON; badge labels should make staleness obvious.

## Additional Repo-Specific Notes
- **Driver script**: `ci/hpc/run_hpc_ci.sh` now exists with `--tests-only`, `--benchmarks-only`, and `--no-codecov` switches. It emits `hpc-ci/hpc-results.json` (schema_version=1) with top-level keys: `commit`, `branch`, `generated_at`, `hpc{cluster,job_id,node,julia_version}`, `coverage{flag,coverage_percent,total_lines,covered_lines,report_file,summary_file|null}`, `tests{status,duration_seconds,log_file}`, `benchmarks{status,duration_seconds,tutorials_log,benchmarks_log}`. Coverage uploads use flag `hpc-extended` by default.
- Keep Julia dependencies centralized in `src/Mycelia.jl`; avoid new `using` statements in leaf files. Any Julia helper for `run_hpc_ci.sh` should live under `src/` or `ci/hpc/` and be included via `Mycelia.jl` if needed.
- Use existing test harnesses: `ci/run_extended_tests.jl` for broad coverage; `run_extended_tests.jl tutorials` / `run_extended_tests.jl benchmarks` are available patterns.
- Prefer `joinpath` for paths and deterministic seeds (per testing guidelines). Reserve HPC walltime in scripts conservatively to avoid scheduler preemption.
- For README badges, follow existing badge style and place new badges near current Codecov badge; document the meaning of each flag in `docs/` or `planning-docs/TOOL_WRAPPER_STATUS.md` if cross-referenced.

## References
- [1]: Jacamar CI Documentation
- [2]: FZ Jülich Apps (Jacamar runner usage examples)
