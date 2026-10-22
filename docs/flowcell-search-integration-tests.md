# Project listing and search integration checks (#163)

Run these with the candidate test image against its normal configuration and
manager directory. Use a read-only manager/config mount for isolated command
checks. No processing daemon, workflow rerun, cleanup or notification is required
for this presentation change. Retain existing mixed-project/correction history;
do not reset the test manager merely to test listing. Substitute known test
identifiers for these synthetic examples.

```console
fm --help
flowcell-manager --help
fm list --help
fm search --help
fm list
fm search GCF-2099-001
fm search 990101_SYNTHETIC_0001_TEST
fm search synthetic
fm search GCF-2099 --status failed
fm search GCF-2099 --stage analysis
flowcell-manager search GCF-2099-001
```

Confirm one project per flowcell row. A project query excludes unrelated projects
on that run; a run query includes its projects subject to row filters. The same
project on two runs appears twice. Both executables must produce identical rows.
Inspect completed, archived, failed and pending history, including a scoped
correction: current analysis/delivery and generation must remain distinct from
last successful delivery. Completed A must stay completed while B is failed or
running, and `--status completed` must still find A despite a failed parent.
Check priority/source against frozen scheduling records, without reading mutable
SampleSheets. See [the column and mapping contract](project-listing.md).

For inventory-only history, verify deduplication by project/normalized run, missing
outputs, per-project archive markers and the `legacy` selectors. If JSON owns a
run, no filter or search may resurrect its stale inventory projects. Older state
must display unavailable detail honestly; an unidentified queued/discovered run
must remain visible with its placeholder.

Compare manager files/directories before and after the commands. Their contents
must be unchanged, including attempts, progress, completions and notifications.
A nonexistent manager directory must remain nonexistent. Read-only mounts should
work. No mail or output changes are expected.

```console
fm search no-such-flowcell-163
echo $?
fm search '   '
echo $?
```

The first prints `No matching flowcells.` and exits 0. The second reports an
invalid query and exits 2. Invalid state reports an error and exits 1.

Automated synthetic coverage in `tests/test_flowcell_search.py` and
`tests/test_project_listing.py` exercises these cases. Installed-command tests use
isolated configuration and mocked SMTP, including subprocesses; the development
gate additionally verifies editable/wheel help outside the source directory.
Record the candidate source/base, build result and actual server command checks
in the PR. Keep operational data and raw transcripts in private evidence storage.
