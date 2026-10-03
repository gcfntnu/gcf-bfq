# Authoritative library-preparation configuration (#123)

BFQ always captures `/opt/gcf-workflows/libprep.config`. Setting
`BFQ_LIBPREP_CONFIG` has no effect and produces a retirement warning. Missing,
unreadable, empty or malformed configuration fails explicitly. The workflow
checkout is copied from its working tree, so local uncommitted edits are included;
the captured configuration bytes then replace the copied `libprep.config`.

## Execution and restart policy

1. Under the existing execution lease, new demultiplexing and analysis attempts
   capture the authoritative configuration once. Basic configuration validation
   happens here, before BCL conversion for new runs.
2. Before launching analysis, shared gcf-tools code selects the kit using actual
   non-index read lengths in output `Stats/Stats.json`. All lanes must agree.
   Kit matching ignores letter case and outer whitespace. The matching ` SE` or
   ` PE` entry wins over an unsuffixed entry; explicit suffixes and layout settings
   must agree with the read geometry. Missing geometry or unsupported kits fail
   within the analysis boundary before configmaker/Snakemake is launched.
3. Every project receives the same captured bytes. Configmaker receives that
   explicit path plus SHA-256, entry and geometry assertions. If the supplied
   snapshot or geometry changes, it fails rather than generating a different
   analysis. Nested kit settings, including `filter.trim`, are retained alongside
   configmaker's project options such as `filter.subsample_fastq`.
4. BFQ keeps the capture in memory. Each project's workflow copy contains the
   exact `libprep.config`, and its generated `config.yaml` contains the
   `libprep_selection` diagnostics (source, SHA-256, kit, entry, geometry and
   workflow). The original authoritative source and selection are also logged.
   There are no separate libprep configuration or manifest files in flowcell
   output. Successful workdir retention is described in the
   [operational README](../README.md#retained-analysis-workdir-snapshots).
5. On successful analysis, BFQ records only the selected workflow name in the
   existing `stages.analysis.metadata.workflow` state field. No state-schema
   version change is needed. Reporting/finalization retries use this name without
   loading `libprep.config` or reconstructing a kit selection. They therefore work
   even if the authoritative configuration has subsequently changed or disappeared.
6. `rerun --from analysis` captures the current authoritative configuration again.
   Existing restart invalidation clears analysis metadata for demultiplexing or
   analysis reruns, and preserves it for reporting/finalization reruns. The new
   workflow is recorded only once the new analysis succeeds.

Older state without a recorded workflow can recover the name from the actual
project `config.yaml` files under `TMPDIR`, then record it in state. Every
identified project's configuration must be readable and agree on the workflow.
Missing, malformed or conflicting project configurations produce an actionable
error: restore the original configs or restart from analysis. BFQ never uses the
current `/opt` configuration to guess how an earlier analysis ran.

The first revision of this PR wrote `bfq-libprep.config` / `bfq-libprep.json` into
flowcell output. Those files are no longer written or read; any files left by an
earlier test build can be removed. They are not required for recovery.

Unknown kits never fall back to an unrelated workflow. For intentional generic
QC, choose a configured entry whose `workflow` is `default`, such as
`Libprep,Custom` (selecting `Custom SE` or `Custom PE`). Correct misspellings in
the effective output SampleSheet or add the appropriate authoritative kit entry,
then queue the normal analysis restart.

The shared `configmaker.libprep` API is used alongside the implemented
[metadata preflight](input-preflight.md) (#121 / gcf-tools#56). Preflight validates
the effective SampleSheet/workbook pair before processing; it does not silently
correct curated metadata or validate every workflow-specific scientific option.

## Build and dependency order

BFQ requires `gcf-tools>=0.3.0` and metadata validation API 1, including the shared
libprep API. BFQ's interpreter and `/opt/conda/bin/configmaker.py` must both use
compatible tools. Production builds select gcf-tools `master`; local development
and CI use the immutable baseline in [requirements-dev.txt](../requirements-dev.txt).
The earlier shared-libprep PR #57 / package 0.2 established this API, but is no
longer a sufficient dependency baseline for current BFQ.

Use the existing [build instructions](../README.md) with your chosen test tag.
For a coordinated tools change, select its actual branch explicitly with `-t`;
do not reuse the historical `bfq-123-shared-libprep-config` branch as a default.

Keep the image's normal `BFQ_ENV=test` setting. Within the image, check:

```bash
/opt/conda/bin/python -c 'from configmaker.libprep import LibprepConfig; from importlib.metadata import version; print(version("gcf-tools"))'
/opt/conda/bin/configmaker.py --help
```

The help must include `--libprep-config`, `--libprep-sha256`, `--libprep-entry`,
and `--expected-read-geometry`.

## Server smoke test

Use a disposable test flowcell, for example `260918_MN00686_0026_A000HCMFHF`, with
the appropriate effective SampleSheet/submission form. Preserve your original
configuration before editing it.

1. In `/opt/gcf-workflows/libprep.config`, change a clearly identifiable parameter
   in the kit's actual SE/PE entry, for example its existing fastp options. Leave
   that edit uncommitted. Optionally set `BFQ_LIBPREP_CONFIG` to a deliberately
   different valid file; expect a warning that it is ignored.
2. Run a fresh test flowcell and inspect the logged source, SHA-256, entry and
   workflow. Check each analysis project's `src/gcf-workflows/libprep.config`
   and generated `config.yaml` (including `libprep_selection`).
   The chosen fastp parameter must appear in `filter.trim.fastp.params`; the
   copied config hash must match the logged hash. Check the actual Snakemake
   fastp command as well, since workflow-level defaults may also affect it.
   After analysis succeeds, check `stages.analysis.metadata.workflow` in the
   flowcell state; no separate `bfq-libprep.*` files should be created.
3. Change the parameter again, then queue a direct analysis rerun:

   ```bash
   flowcell-manager rerun 260918_MN00686_0026_A000HCMFHF --from analysis
   ```

   Confirm the new hash and parameter reach the regenerated configuration and
   reports. Also check that an existing project's FASTQ directory is reused on
   this rerun (the `--skip-create-fastq-dir` check now uses the project directory).
4. A reporting/finalization retry must use the completed analysis's recorded
   workflow, even if `/opt/gcf-workflows/libprep.config` has since changed or is
   unavailable. It must not read or create a separate retained libprep file.
5. On the disposable run, exercise an unknown kit and malformed/missing config.
   Expect an actionable error report and no affected analysis subprocess launch.
   Restore valid inputs/configuration and retry through normal manager commands.

Local verification covers the BFQ execution path and real configmaker API/CLI
with temporary fixtures. Actual BCL conversion, Snakemake/container execution,
server mounts and notification delivery still require the server test above.

## Explicit analysis resume

`fm rerun --from analysis --resume` preserves the retained execution context.
Shared input validation still applies; curated identifiers/wells must agree with
the retained configuration. It does not regenerate configmaker output or read
installed libprep settings. See [analysis resume](analysis-resume.md) for the
binding checks, operator correction path and server verification.
