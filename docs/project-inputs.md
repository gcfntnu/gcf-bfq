# Project-paired flowcells and scoped corrections

A flowcell can contain independent projects using different library preparations,
organisms, read layouts and workflows. Supply one SampleSheet and one submission
form per project. For example:

```text
SampleSheet_GCF-2099-001.csv
GCF-2099-001_Sample-Submission-Form.xlsx
SampleSheet_GCF-2099-002.csv
GCF-2099-002_Sample-Submission-Form.xlsx
```

SampleSheet names start with `SampleSheet` and end in `.csv`; workbook names
contain `Sample-Submission-Form` and end in `.xlsx`. Each filename must contain
exactly one `GCF-YYYY-NNN` project key. File order, timestamps and overlapping
sample names do not establish pairing. Each sheet's `Sample_Project` must agree
with its filename. Historical single-pair inputs keep their existing filenames
and may still describe multiple projects.

Exact duplicate SampleSheets are ignored deterministically before pairing,
preferring a project-qualified filename over a generic copy. Their names and
SHA-256 hashes are recorded. Identical bytes under filenames naming different
projects are an error. Missing partners, duplicate forms, conflicting sheets
and unequal project-key sets must be corrected explicitly.

An explicit whole-run refresh can replace a legacy canonical pair with project
pairs. Obsolete canonical input bytes are retained under
`provenance/input-replacements/` so they do not become competing candidates.

Every effective pair is validated independently by `gcf-tools`, followed by
whole-flowcell lane, index, mismatch-policy, sample/output identity and physical
read-geometry checks. A different project's workbook cannot supply missing
metadata. All projects must pass even during a correction selecting one project.
The existing extra-sample, ASCII identifier and pooled-library rules still apply.

Compatible non-10x projects share a structured BCL Convert pass. Uniform settings
use a v1 sheet with legacy statistics; per-sample settings use a v2 sheet with
native conversion reports. Derived v2 sheets omit the optional `Sample_Name`
column so FASTQ filenames continue to use `Sample_ID`; original names remain in
the source sheets and row provenance. Enabling name-based output would change
the output layout ([Illumina sample-sheet reference](https://support-docs.illumina.com/SW/dragen_v42/Content/SW/DRAGEN/SampleSheet.htm)).
Physical read geometry remains available from the
frozen `RunInfo.xml` when native reports omit `Stats.json`. The supported
per-sample settings are `OverrideCycles`, `BarcodeMismatchesIndex1`,
`BarcodeMismatchesIndex2`, `AdapterRead1`, `AdapterRead2`, `AdapterBehavior` and
`AdapterStringency`; this requires BCL Convert 4.1 or newer (the workflow image
currently specifies 4.2.4). Other setting conflicts fail with an actionable
diagnosis before cleanup. An incompatible backend combination also fails;
multi-pair 10x backend/index normalization belongs to the linked #156 work.
`FORCE_BCL2FASTQ` remains available for legacy single-pair runs and is explicitly
rejected for multi-pair inputs. There is no silent backend fallback or
scientific-setting substitution.

The authoritative library-prep configuration is captured once per attempt.
Each project's kit resolves against its own effective biological reads; physical
sequencing geometry remains separately recorded. Structured conversion provenance
retains derived sheets, source row mappings, hashes, effective settings and tool
identity. Analysis configuration, delivery inputs and retained snapshots identify
their own project and execution attempt.

Optional SampleSheet `Priority` controls deterministic ordering within analysis
and archiving, with all selected analyses finishing before archiving starts.
See [project scheduling](project-scheduling.md) for decimal values, defaults,
manager previews and the frozen execution plan.

## Correct one project

`fm` and `flowcell-manager` remain the same installed interface. Preview before
queueing a correction, repeating `--project` to select multiple projects:

```bash
fm rerun RUN_ID --from analysis --project GCF-2099-002 --refresh-inputs --dry-run
fm rerun RUN_ID --from analysis --project GCF-2099-002 --refresh-inputs
fm rerun RUN_ID --from demultiplexing --project GCF-2099-002 --refresh-inputs --dry-run
```

Omitting `--project` retains whole-run rerun behavior. Curated output-side inputs
win independently per project. `--refresh-inputs` replaces only the selected
projects' inputs in a scoped correction. The preview identifies selected and
retained projects, proposed input replacements, restart stages and cleanup scope.
`--force` skips confirmation and does not bypass validation or execution leases.

Selected corrections use frozen input copies and attempt-owned working outputs.
An analysis correction preserves FASTQs and valid FASTQ manifests. A conversion
correction can repeat the combined pass into staging; only selected projects'
replacement FASTQs and dependent products are published. Regenerated unselected
products are not installed. Unselected FASTQs, timestamps, checksums, analysis,
reports, archives and successful snapshots remain intact.

The operator decides which projects require correction and must widen the scope
when appropriate. BFQ does not require read-assignment-equivalence proofs or
automatically diagnose sample mix-ups. A failed correction retains the prior
successful delivery and provenance. Publication uses the existing state store
and execution lease, with a journal for recovery around the processing completion
commit. State distinguishes historical project completion from the current
attempt. Notifications remain independent of processing success.

Flowcell-wide statistics describe a conversion attempt, rather than a sum of
repeated passes. A scoped correction retains its attempt statistics/provenance
without overwriting those associated with retained projects. Inspect the state,
attempt records and correction provenance when identifying the source of a report.

## Verification boundaries

The local Python gate exercises synthetic pairing, validation, planning, scoped
publication and failure recovery, plus legacy regression tests. External
processing and mail are mocked in those checks. They do not establish scientific
output correctness. Facility integration must run the mixed-workflow fixture,
inspect each project's reports/archives/snapshots and checksums, then exercise
both scoped analysis and scoped conversion corrections while verifying the
retained project's bytes and timestamps. Candidate image/source identities and
actual performed integration checks belong in the implementing PR; raw facility
data and operational records remain private.
