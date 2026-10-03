# Resume an existing analysis (#143)

Ordinary `fm rerun RUN_ID --from analysis` still deletes downstream products and
rebuilds configuration/workflows from curated inputs and the installed workflow.
For recovery that should reuse completed scientific jobs, explicitly select:

```bash
fm rerun RUN_ID --from analysis --resume --dry-run
fm rerun RUN_ID --from analysis --resume --reason "Continue after local repair"
```

`fm` and `flowcell-manager` are the same entry point. `--force` only skips manager
confirmation; it never forces scientific recomputation. Resume requires an
explicit analysis boundary and rejects `--refresh-inputs` and index-correction
options. It does not initialize legacy inventory-only runs.

## What is preserved and what runs

The manager preview is read-only, shows every retained workdir, configuration and
workflow, and describes the state transition. It is **not a Snakemake DAG dry
run**. Queueing removes no outputs: `data/`, `.snakemake/`, FASTQ links, config,
Snakefile, PEP, workflow and libprep files, existing delivery products and previous
successful snapshots all remain. Demultiplexing, valid FASTQ manifests and early
sequencing QC remain valid. Analysis becomes queued; reporting and finalization
become pending. Downstream notification intents are superseded normally.

Execution uses those workdirs without running configmaker, copying installed
workflows, or selecting today's `/opt/gcf-workflows/libprep.config`. It invokes
the usual `multiqc_report` target with normal dependency/provenance triggers and
`--rerun-incomplete`. A complete project may legitimately be a
no-op. In a multi-project run, **every** project must have a usable owned workdir;
there is no destructive fallback for an absent project.

After successful Snakemake execution, BFQ refreshes report and metadata copies
and replaces the project's output-side QC copy (including removal of stale
members). Reporting and archive creation/checksumming then run normally, followed
by state completion, retained snapshot publication and notification handling.
Queueing and failed scientific execution retain previous delivery files; those
files are not evidence of current completion. Publication across all projects is
not a transaction: an earlier successful project's report copy can be refreshed
before a later project fails. The previous successful retained snapshots remain
protected until normal successful finalization.

Snakemake retains its **normal failed-job output cleanup**. BFQ does not use
`--keep-incomplete`: the real tiny-DAG test exposed that Snakemake 9.7.1 can mark
a failed job's partial output as complete with that flag and return a false no-op
on an unchanged retry. Completed reusable branches remain untouched; outputs
belonging to the failed job may be removed by Snakemake. No BFQ implementation
edits `.snakemake` metadata or substitutes a blanket force/recompute option.
The production Snakemake pin is unchanged.

## Validation and local repairs

This is reuse of an operator-reviewed execution context, not regeneration from
metadata. BFQ still runs the shared gcf-tools input validator before queueing and
execution. Both canonical output-side files must exist and validate. An invalid
original workbook does **not** get a resume exemption.

For the `ÅLE01` incident, correct the output-side Cell Multiplexing IDs as well as
the four retained config `wells` keys and embedded `Sample_ID` values. Keep `Wells`
assignments consistent. Do not regenerate configmaker outputs just to apply that
repair. Record the reason when queueing. `fm validate RUN_ID` checks the workbook
and SampleSheet; the resume preview additionally checks the retained context.

Resume checks the BFQ ownership marker against the full run ID/project and any
recorded ownership token, retained config project/sample/flowcell identities,
PEP sample identities and required files, current parsed Cell Multiplexing
ID/well assignments when present, and FASTQ symlinks covering this output
run/project. All projects must agree on the retained workflow. Unmarked older
workdirs cannot prove ownership and are refused; do not fabricate ownership
markers to bypass this. Restore an owned workdir or use an ordinary rerun.

The preview/queue binds source input SHA-256 hashes, config/Snakefile/PEP/workflow
file hashes (including local edits, excluding `.git` and Python bytecode caches),
workdir ownership and FASTQ link targets/sizes/mtimes. BFQ rechecks the binding
under its execution lease at queue and execution. If it changed, processing fails
with an instruction to inspect and requeue; no config regeneration occurs.
Make intentional repairs **before** queueing. Do not edit or run the workdir
concurrently with BFQ. The per-flowcell lease does not protect a separate manual
Snakemake process or another flowcell sharing the same project/date directory.

The binding is stored in `restart_request.analysis_resume` and copied to
`attempts[].analysis_resume`, so it survives daemon restart and remains in attempt
history. The attempt also records the command, resolved Snakemake executable and
BFQ interpreter at execution. Installed runtime versions remain recorded as runtime information;
the resume binding and successful workdir snapshot describe the actual retained
analysis source, including local modifications. Existing timing machinery
records only elapsed work in the new attempt; it does not invent durations for
reused jobs or turn a previously failed analysis duration into successful time.

BFQ does not validate every scientific option or infer dependencies omitted from
the DAG. FASTQ identity checks do not reread large FASTQs to prove content
identity. Use an ordinary invalidating rerun when the scientific change requires
recomputation Snakemake cannot infer. Metadata corrections beyond the checked
identities/well assignments still require operator review of retained config/PEP.

## Local evidence and remaining server checks

`python3.11 scripts/dev.py check all` includes:

- Unit coverage for read-only preview, no-delete queueing, restart persistence,
  ownership/config/PEP/link failures, input binding, invalid options, active leases,
  multiple projects and report refresh without configmaker.
- The guarded installed-command operational scenario: preserved repairs and
  expensive-result/provenance bytes and mtimes, failed resume retaining the old
  snapshot, then a doubled Snakemake no-op completing normal BFQ delivery despite
  a doubled SMTP failure. No scientific commands run in this profile.
- A separate real **tiny local** Snakemake 9.7.1 DAG in the fast profile. Its two
  shell jobs write short text files; it tests a downstream failure/recovery,
  preserved expensive-branch bytes/mtime, no-op and ordinary parameter-triggered
  reruns. Its cache and outputs are invocation-local. This is the only real
  workflow execution in the routine checks; it uses no containers or genomics.

Before merge, use disposable operational data in the BFQ test image:

1. Retain a partially completed analysis and its valid FASTQ manifest. Back up and
   apply a bounded repair consistently to curated metadata/config as needed. Save
   checksums/mtimes of a completed expensive branch and the previous snapshot.
2. Preview both command aliases; confirm no deletion or file/state changes. Run
   an actual Snakemake dry run separately in each workdir to inspect the DAG.
3. Queue resume, restart BFQ before it executes, and confirm the mode survives.
   Confirm configmaker/workflow copying never occurs, repaired bytes remain, the
   expensive branch is reused and only the DAG-selected jobs execute.
4. Repeat with a deliberately changed installed workflow/libprep; retained source
   must still be used. Change a DAG-visible parameter locally before requeueing
   and verify the appropriate dependent jobs run.
5. Exercise a multi-project run (one already complete), missing/mismatched workdir,
   active lease, and a change after queueing. Unsafe cases must stop without
   deleting outputs or refreshing configuration.
6. After manual completion, exercise a no-op resume through actual report copying,
   MultiQC, 7-Zip, archive checksums, inventory/state and snapshot publication.
   Inspect removed stale report members, archived content and retained local edits.
7. Fail a resumed job, verify preserved work/old snapshot, repair and requeue.
   Check early QC/manifests and notification supersession. Actual email delivery
   requires a separately arranged test; local tests never send mail.
8. Verify an ordinary analysis rerun still resets work and selects installed code.

Merge, production promotion and deployment remain separate maintainer decisions.
