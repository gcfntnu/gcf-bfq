# SampleSheet and submission-form preflight

BFQ validates the effective `SampleSheet.csv` and `Sample-Submission-Form.xlsx`
before demultiplexing and before analysis, including an analysis-only restart and
restored FASTQs. Domain rules and workbook parsing come from the same gcf-tools
validator used by standalone configmaker.

Every SampleSheet sample must exist in the **effective merged** submission form.
Extra submission-form samples are allowed. Repeated lane rows for the same sample
are allowed when their assignments agree. Unreadable workbooks, unusable headers,
ambiguous sample mappings, and samples lost during the customer/lab worksheet
merge produce actionable errors. Informational findings and warnings alone do
not block a run. Identifiers and input contents are never repaired automatically.

BFQ also requires a known sequencing instrument before demultiplexing or analysis.
It resolves the second underscore-separated field of the run ID through
configmaker's `SEQUENCERS` registry in gcf-tools. A missing instrument field or
unregistered instrument produces a `run.sequencer` preflight error naming the run
and, when available, the instrument ID. Check the run identity and the installed
gcf-tools mapping before retrying; a new instrument must be registered in
gcf-tools. BFQ does not guess a model or continue with workflow defaults.

This check applies to `validate`, initialization, processing reruns (including
project corrections and `--resume`), and execution. `--force` does not bypass it.
Rerun refusal preserves existing results, FASTQs, checksums and successful
snapshots. Delivery-only reporting/finalization retries retain their existing
behavior. Historical notification summaries display `Unknown` when lookup fails;
notification replay does not require a currently registered instrument.

The non-processing `get_sequencer()` lookup returns `""` and logs a warning for an
unresolved instrument, matching configmaker's unknown-value convention. BFQ's
processing guard rejects that condition before configmaker or Snakemake starts,
so diagnostic prose is never passed as `--machine`. Known machine names are
passed unchanged. Instrument mappings and scientific optical-duplicate settings
remain owned by gcf-tools and gcf-workflows respectively.

## Check inputs without changing a run

```sh
flowcell-manager validate RUN_ID
```

This works before discovery, demultiplexing, or analysis outputs exist. It prints
the selected paths and the complete findings; exit status is zero on success and
nonzero on error. It does not create state directories, acquire execution leases,
recover interrupted state, queue processing, copy inputs, or remove results.

Each canonical output-side input is authoritative, including malformed files. If
only one output-side input exists, BFQ preserves it and selects the missing partner
from the instrument. It no longer replaces both files when only one is missing.
Legacy alternative filenames are normalized by copying their exact bytes during
execution or rerun preparation. On the instrument, SampleSheets with BFQ
`[CustomOptions]` retain their discovery precedence; ordinary instrument sheets
without those options do not opt a new run into automatic BFQ processing.

To preview an explicit instrument refresh without copying anything:

```sh
flowcell-manager validate RUN_ID --refresh-inputs
```

This selects both instrument inputs, matching the pair a rerun with
`--refresh-inputs` will install. The command can validate a running run without
modifying its state, but the result describes the bytes read during the check,
not a promise about concurrent edits. The daemon always validates again.

## Correct and retry

1. Read the error report or run `flowcell-manager validate RUN_ID`. Findings name
   the affected file, worksheet, row, column, and sample where available.
2. Correct the effective output-side file(s), keeping intended IDs and metadata.
   If the instrument copies are the intended replacement, correct those instead
   and use the explicit refresh commands above and below.
3. Validate the corrected pair, then queue the existing restart boundary:

   ```sh
   flowcell-manager rerun RUN_ID --from demultiplexing
   # For a failure at the analysis boundary with retained FASTQs:
   flowcell-manager rerun RUN_ID --from analysis
   # Explicitly replace both effective inputs from the instrument:
   flowcell-manager rerun RUN_ID --from analysis --refresh-inputs
   ```

Choose the boundary at which the run failed, or an earlier applicable boundary.
Changing demultiplexing assignments (sample IDs, projects, indexes, or lane
assignments) requires restarting demultiplexing. Correcting metadata while keeping
the FASTQ assignments intact can use analysis. Check the command's cleanup preview:
an analysis rerun removes downstream reports, analysis summaries, archives, and the
analysis work directory while retaining FASTQs and inputs.

Rerun and explicit initialization validate before cleanup. Invalid inputs leave
existing downstream results in place; correcting them needs no state-file edits.
The command checks again under the execution lease, installs effective inputs,
and verifies the installed pair before removing results. A copy or cleanup failure
can leave the run in `preparing`, which requires an explicit operator retry.
`--dry-run` performs the input check but makes no changes. `--force` bypasses the
confirmation prompt, not validation or the execution lease.

## State, reports, and notifications

Input preflight is part of the existing `demultiplexing` or `analysis` boundary;
there is no additional public restart stage. The daemon performs it under the
execution lease before expensive processing and before FASTQ manifest generation
on direct analysis entry.

The complete structured report includes SHA256 hashes of the parsed inputs, the
validator version, findings, planned sample/project/Sample_Group summaries, and
the check time and restart boundary. BFQ stores it in:

- `<reportDir>/<RUN_ID>.input-preflight.json` for the latest execution check;
- `stages.<boundary>.metadata.input_preflight` in `flowcell-manager show RUN_ID`;
- the corresponding attempt's `input_preflight` history.

Every execution revalidates current bytes. A previous successful report never
authorizes changed files or a different validator version. Existing error-report
and email routing is retained. A preflight failure identifies the check explicitly
and includes the complete diagnostics. Planned input samples remain distinct from
the samples configmaker actually discovers in FASTQs during analysis.

## Server integration checks

Build both BFQ and its configmaker environment with the coordinated gcf-tools
version before testing. Use a disposable run/output copy and the test notification
route.

1. Put a Sample_ID in the SampleSheet that is absent from the merged submission
   form. `validate` must fail; start BFQ and verify that no demultiplexer process
   launches. State should be failed at demultiplexing, with a hash-bound report
   and a useful error email.
2. Correct the output-side workbook. `validate` must pass; rerun from
   demultiplexing and verify normal processing. Change the instrument copies to
   disagree first if needed to confirm that the curated output pair wins.
3. With completed outputs present, make the effective pair invalid and request an
   analysis rerun. It must fail before deleting QC reports, archives, or workdirs.
   Correct the pair and rerun; verify that the daemon checks again before analysis.
4. Repeat the preceding check using invalid instrument inputs and
   `--refresh-inputs`; it must preserve the existing output pair/results on failure.
   A successful refresh must use the instrument pair for execution.
5. Restore a recognized FASTQ tree containing the input pair, then test an invalid
   pair and a corrected retry from analysis. No demultiplexer should launch.
6. Check a valid subset with additional submission-form samples and legitimate
   repeated lane definitions. Both manual and automatic validation must pass and
   retain useful sample/project/Sample_Group summaries.
7. Test malformed workbook/header data and a customer sample missing from a
   nonempty lab sheet. Reports/emails must remain readable. Repeat an unchanged
   failure to check duplicate suppression, then alter the offending input bytes
   and verify the changed failure gets a new notification identity.
8. Compare planned samples in the preflight report with configmaker's actual
   FASTQ-discovered analysis summary; missing FASTQs must not be reported as
   successfully analyzed merely because their inputs passed validation.
9. On a disposable fixture with an unregistered synthetic instrument ID, verify
   that `validate` and forced demultiplexing/analysis reruns fail before cleanup.
   Execution must record `run.sequencer` in the failed stage's preflight report
   without launching conversion, configmaker or Snakemake. Repeat with a missing
   instrument field and confirm an actionable diagnostic. On a registered
   fixture, verify the expected machine name in the generated configmaker config
   and normal workflow startup. Use the test notification route to confirm that
   historical messages with an unknown instrument still display `Unknown`.

## Explicit analysis resume

`fm rerun --from analysis --resume` preserves the retained execution context.
Shared input validation still applies; curated identifiers/wells must agree with
the retained configuration. It does not regenerate configmaker output or read
installed libprep settings. See [analysis resume](analysis-resume.md) for the
binding checks, operator correction path and server verification.
