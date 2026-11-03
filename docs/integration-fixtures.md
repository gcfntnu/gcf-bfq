# Operational integration fixtures

This guide describes the BFQ contract for a maintained server fixture setup.
The [development checks](development.md) exercise synthetic data and mocked
external boundaries. Scientific processing, container compatibility and actual
notification delivery require a separately arranged server integration run.

The site owns fixture selection, input storage, launchers, resource overrides
and credentials in its private operations repository/runbook. Those details are
not prerequisites for understanding or developing BFQ. BFQ itself provides the
normal discovery, initialization, state and restart interfaces used below.
Site-specific fixture selection helpers are not installed BFQ commands.

## Define coverage and starting points

Record the workflow, input kind, starting stage and expected acceptance checks
for each fixture. A routine profile can select a representative case per
supported workflow, plus distinct orchestration cases such as multiple projects
and materially different demultiplexing routes. Keep additional cases available
for changes that affect them. Agree the required scope before a release; a
routine subset does not establish coverage of every maintained case.

| Input | Start | What it establishes |
| --- | --- | --- |
| Raw instrument data with valid run/project metadata | `demultiplexing` | Conversion, downstream analysis, reporting and finalization |
| Restored FASTQs with valid SampleSheet/workbook inputs | `analysis` | Analysis and downstream behavior; no conversion validation |

Use the supported initializer after materializing inputs into an isolated test
installation. These examples use a placeholder run ID; substitute the selected
test run and its configured input path:

```console
fm initialize /instruments/RUN_ID --from demultiplexing --dry-run
fm initialize /instruments/RUN_ID --from demultiplexing
```

For restored FASTQs, place the validated project FASTQs and curated metadata in
the configured output tree and initialize from `analysis`. Follow the
[input authority and initialization rules](../README.md#flowcell-manager).
Use `fm rerun`, with a deliberate restart boundary, for an existing state-backed
run. Repeated setup must preserve selected, running and completed states; it
must not silently queue completed work again.

## Isolate mutable state

Each independent test session needs its own discovery/import roots, output tree,
manager, analysis workdirs and temporary paths. Keep original fixture sources
read-only and materialize independent mutable copies when required. Retain prior
sessions and failed outputs for diagnosis. Do not reset shared reference or
container-image caches as part of session selection. Coordinate cache population
and resource use when more than one installation shares them.

An empty test selection requires **both** empty discovery/import roots and no
queued state or eligible restored output. Clearing only the manager can make
existing instrument directories eligible again. Project holds or fabricated
completion records are not a substitute for an empty selection.

Stop the relevant daemon and check execution leases before switching session
paths. Use supported manager commands for lifecycle changes; never edit state
JSON or delete successful outputs to obtain a passing result.

## Reduce fixtures without losing their purpose

Preserve paired-read order and identity, metadata relationships and valid MD5
manifests. Record original/derived checksums and a deterministic reconstruction
recipe. Derivation should depend on stable source inputs rather than a previous
test's output files.

For WGS, a separately named assembly containing intact reference sequences can
reduce indexing/alignment cost while retaining sequence names and coordinates.
Retain corresponding read pairs and verify meaningful mapped coverage. For
classification, select a smaller database through the normal workflow
configuration and retain enough reads for useful classification and abundance
results. Small fixtures test integration; they are not accuracy benchmarks for
the third-party scientific tools.

Scientific defaults, reference getters and library-prep configuration remain
owned by `gcf-workflows`; shared metadata selection remains owned by `gcf-tools`.
Apply test-specific overrides only to the test installation, preserving the
normal configuration-authority rules. Distinguish resource preparation from
warm execution when comparing runtimes.

## Acceptance and evidence

Wait for terminal processing **and finalization** for every selected run. Check:

- Expected projects, samples, starting stages and actual selected workflows.
- FASTQ integrity and MD5 manifests; conversion-specific read oracles where
  available for raw-data fixtures.
- Meaningful workflow outputs, including route-specific numerical checks.
- Reports, delivery archives and retained successful analysis snapshots.
- Notification coverage and SMTP outcomes; distinguish SMTP acceptance from
  independently observed mailbox receipt.
- Empty selection, single selection, repeated setup and completed-state restart
  behavior when the fixture setup itself changes.

Record fixture/setup revision, actual source/dependency identity, selected
scope, attempts, timings and limitations. Keep sensitive evidence in the site's
private store. Local tests must mock SMTP, including subprocesses; real server
notifications belong only to the arranged integration run. Passing tests does
not itself authorize a merge, production promotion or deployment.
