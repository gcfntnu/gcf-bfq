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
4. The flowcell output retains `bfq-libprep.config` and `bfq-libprep.json` with
   source, hash, input kit, entry, geometry and workflow. The same selection is
   visible in BFQ logs and `bcl2fastq.ini`; each project `config.yaml` also records
   the snapshot it consumed under `libprep_selection`.
5. `rerun --from analysis` captures the current authoritative configuration again.
   Reporting/finalization restarts restore the prior analysis selection, even if
   the authoritative source is subsequently edited or removed. Corrupt or partial
   retained snapshot pairs fail with a repair/restart instruction.

Existing results without either snapshot file are bootstrapped from the current
authoritative configuration and existing Stats.json with an explicit warning.
Check that the logged workflow matches those legacy results; use an analysis
restart if configuration changes should apply. There is no retained history of
every configuration attempt, and workflow code is still copied per project.
A complete workflow/provenance snapshot remains separate future work.

Unknown kits never fall back to an unrelated workflow. For intentional generic
QC, choose a configured entry whose `workflow` is `default`, such as
`Libprep,Custom` (selecting `Custom SE` or `Custom PE`). Correct misspellings in
the effective output SampleSheet or add the appropriate authoritative kit entry,
then queue the normal analysis restart.

The shared `configmaker.libprep` API is portable and can be called by the future
preflight work in #121 / gcf-tools#56. This change does not implement workbook
validation, metadata correction, or broader workflow-specific parameter checks.

## Build and dependency order

Companion dependency: [gcf-tools PR #57](https://github.com/gcfntnu/gcf-tools/pull/57),
branch `bfq-123-shared-libprep-config`, package version 0.2. BFQ now requires
`gcf-tools>=0.2`. Both BFQ and `/opt/conda/bin/configmaker.py` must use that version.
BFQ CI pins the tested companion commit so it can verify this PR before merging
the dependency. Production promotion must include gcf-tools in `master` before
building BFQ's updated production image.

From this BFQ issue branch, build locally for server testing:

```bash
docker build -f dockerfile-test -t gcfntnu/bfq:issue-123 \
  --build-arg GCF_TOOLS_BRANCH=bfq-123-shared-libprep-config .
```

Alternatively, the existing build-and-push helper accepts the branch via
`bash build-tag-push.sh test YOUR_TEST_TAG -t bfq-123-shared-libprep-config`.
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
   workflow. Check flowcell `bfq-libprep.json` / `bfq-libprep.config`, each analysis
   project's `src/gcf-workflows/libprep.config`, and generated `config.yaml`.
   The chosen fastp parameter must appear in `filter.trim.fastp.params`; the
   copied config hash must match the logged hash. Check the actual Snakemake
   fastp command as well, since workflow-level defaults may also affect it.
3. Change the parameter again, then queue a direct analysis rerun:

   ```bash
   flowcell-manager rerun 260918_MN00686_0026_A000HCMFHF --from analysis
   ```

   Confirm the new hash and parameter reach the regenerated configuration and
   reports. Also check that an existing project's FASTQ directory is reused on
   this rerun (the `--skip-create-fastq-dir` check now uses the project directory).
4. For a run with two projects, edit the source after the capture log appears.
   Both projects must retain the captured configuration. A later analysis rerun
   must pick up the edit. A reporting/finalization retry must retain the completed
   analysis's workflow even if the source now specifies a different workflow.
5. On the disposable run, exercise an unknown kit and malformed/missing config.
   Expect an actionable error report and no affected analysis subprocess launch.
   Restore valid inputs/configuration and retry through normal manager commands.

Local verification covers the BFQ execution path and real configmaker API/CLI
with temporary fixtures. Actual BCL conversion, Snakemake/container execution,
server mounts and notification delivery still require the server test above.
