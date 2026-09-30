# Flowcell search and fm smoke checks (#116)

Run these in the rebuilt test image against its normal `/config/bcl2fastq.ini`
and manager directory. Search and list do not alter run records, queue work or
send notifications. Replace the example identifiers with runs in that inventory.

```console
fm --help
flowcell-manager --help
fm search --help
fm search GCF-2026-043
fm search 260925_NB501038_0281_AHL2T7AFXC
fm search hl2t7afxc
fm search GCF-2026 --status failed
fm search GCF-2026 --stage analysis
flowcell-manager search GCF-2026-043
```

Confirm that project search returns every associated run, that flowcell search
includes all projects on that run, and that both executables produce the same
rows. Look up one completed and one archived run as well. If JSON and inventory
both describe a run, it should appear once using JSON projects/status/stage;
a `completed`, `archived` or `legacy` filter must not resurrect an old inventory
row for a currently failed/queued run. `fm list` remains the unfiltered overview.

Check exit behavior with an identifier known to be absent:

```console
fm search no-such-flowcell-116
echo $?
fm search '   '
echo $?
```

The first prints `No matching flowcells.` and exits 0. The second reports an
invalid query and exits 2. No real rerun, cleanup or notification is needed to
validate this change: both names resolve to the same existing CLI entry point.

Automated coverage in `tests/test_flowcell_search.py` exercises canonical and
legacy records, literal matching, deduplication, filters and installed command
behavior with isolated configuration. The development verification also installs
the built wheel and runs its help/search checks away from the source directory.
