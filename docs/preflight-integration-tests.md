# Coordinated input validation: integration tests

These checks cover gcf-tools issue 56 and BFQ issue 121. Run them in the test
deployment with a disposable copy of a flowcell and its outputs. Do not merge
either PR until the server checks pass. No production inputs need to be edited.

## Build both changes together

Check out BFQ branch `121-input-preflight`, then use the standard build/tag/push
wrapper with its companion tools branch (choose your intended test image tag):

```bash
bash build-tag-push.sh test preflight-test -t 56-shared-input-validation
```

The wrapper builds and pushes `gcfntnu/bfq:preflight-test`; it defaults the workflow
branch to `bfq-dev`. Use the normal test deployment's mounts and configuration.
The build checks that `/opt/conda/bin/python` has validation API 1 and imports BFQ.
The clean CI environment retains a full `pip check`. Image builds do not apply
that audit to all unrelated applications in the shared Conda environment.
The test image sets `BFQ_ENV=test`, which suppresses production error emails.
Inside the container, verify the tools used by BFQ and configmaker:

```bash
/opt/conda/bin/python -c 'from configmaker.validation import VALIDATOR_VERSION, VALIDATION_API_VERSION; print(VALIDATOR_VERSION, VALIDATION_API_VERSION)'
/opt/conda/bin/python /opt/conda/bin/configmaker.py --help
flowcell-manager validate --help
```

Expected: gcf-tools `0.3.0`, API `1`, and configmaker exposes
`--expected-validation-version`. BFQ passes this version to configmaker, so a
different subprocess installation fails explicitly before project initialization.
For installations outside this image, update gcf-tools in both Python environments.

## Manual validation and input selection

Use an existing test run ID, for example:

```bash
RUN_ID=260918_MN00686_0026_A000HCMFHF
flowcell-manager validate "$RUN_ID"
flowcell-manager validate "$RUN_ID" --refresh-inputs
```

The command prints the two selected paths, findings with codes/context, and the
planned sample summary. Check its exit status immediately: `echo $?` is `0` on
valid inputs and `1` on errors. Warnings and extra submission samples pass.
Neither command copies inputs, creates state, queues work, nor cleans outputs.
`--refresh-inputs` here only previews validation of the instrument pair.

Verify these cases before starting the daemon:

| Input case | Expected result |
| --- | --- |
| SampleSheet is a subset of submission metadata | Pass; extra samples are informational |
| Same sample occurs on distinct lanes | Pass with one unique planned sample |
| SampleSheet ID absent from submission metadata | Error naming the ID and file |
| Customer sample absent from a nonempty lab sheet | Error naming the sample and missing worksheet |
| Empty lab worksheet | Pass if the customer sheet supplies the metadata |
| Populated worksheet row without an ID | Error with worksheet and row |
| Invalid workbook bytes or missing required header | Readable parsing/header error, no traceback cascade |
| Leading-zero ID or surrounding whitespace mismatch | Exact comparison; no automatic correction |
| Curated output pair differs from instrument pair | Default uses output; explicit refresh uses instrument |
| Only one curated output input exists | Preserve it; select only the missing partner from the instrument |

## Execution-time blocking and correction

1. On a fresh test flowcell with valid BFQ `[CustomOptions]`, make the submission
   form incompatible (for example, remove one required sample from the lab sheet).
   Let the daemon discover it.
2. Confirm failure is recorded within **demultiplexing**, with **Input preflight**
   in the log/error. There must be no demultiplexer launch or FASTQ production.
3. Inspect `<reportDir>/<RUN_ID>.input-preflight.json` and the stage metadata
   displayed by `flowcell-manager show "$RUN_ID"`. The report includes
   input paths, SHA-256 hashes, validator version, and all contextual findings.
4. Correct the **effective output-side** input shown by validation, then run:

   ```bash
   flowcell-manager validate "$RUN_ID"
   flowcell-manager rerun "$RUN_ID" --from demultiplexing --force
   ```

   Alternatively, correct the instrument pair and explicitly use
   `rerun ... --refresh-inputs`. Wake/start the daemon using your normal process.
5. Confirm the new report contains the corrected hashes and the run proceeds.
   No manual JSON state edits should be necessary.

Also test the gap between preparation and execution: stop the daemon, queue a
valid rerun, then make its effective pair invalid before starting the daemon.
The execution check must still block it. A previously passed report must never
authorize changed bytes.

## Rerun safety, direct analysis and restored FASTQs

On a completed disposable run, preserve representative downstream products
(MultiQC report, archive and analysis work directory) and note their hashes.
Make the effective pair invalid, then request:

```bash
flowcell-manager rerun "$RUN_ID" --from analysis --force
```

Expected: validation failure before downstream cleanup; the existing products
and their contents remain. Repeat with `--from demultiplexing` and with an
invalid instrument pair plus `--refresh-inputs`.

Correct the inputs and queue analysis again. Confirm validation occurs before
FASTQ-manifest checking/configmaker, the demultiplexer is not rerun, and the
analysis uses the corrected metadata. Repeat on an independently restored FASTQ
run recognized by BFQ's normal restore rules. Metadata validation itself does
not require FASTQs, so manual validation also works before restoration.

## Reports and notifications

- Completion mail must show **Planned input samples (current metadata)**,
  unique per-project counts and Sample_Group values/missing counts.
- A separate section reports **Samples discovered in FASTQs at analysis
  initialization** from `configmaker-analysis-<project>.json`. Use a test project
  with a planned sample lacking FASTQs to verify the counts remain distinct.
- Malformed metadata must produce readable findings without breaking mail
  composition. Error mail includes the normal `.error` report and the structured
  JSON result. The JSON must describe the failed bytes even if the file is edited
  before delivery.
- Test actual error email delivery only in an isolated test configuration with
  `error_to` pointing to the intended test mailbox and `BFQ_ENV=production` for
  that test process. No emails are sent by the automated tests. An unchanged
  failure remains suppressed across retries; changed input bytes/findings produce
  a distinct notification identity. Existing routing is unchanged.

## Unexpected compatibility findings

- The first integration image build installed BFQ successfully but then failed
  the newly added global `pip check`: `snakemake-interface-common 1.23.1` requires
  `packaging>=26.1`, while the image has `packaging 25.0`. The Dockerfiles now check
  BFQ imports and the validator API instead. This removes the unrelated build
  gate; it does not repair the Snakemake dependency mismatch. That belongs to the
  base-image dependency maintenance and should be considered if Snakemake fails.
- Old BFQ could overwrite a curated partner when the other input was absent.
  Each output-side input now takes precedence independently. A malformed curated
  canonical file fails visibly instead of triggering an instrument fallback.
- A historical submission fixture contains a free-text note inside the sample
  table without an ID. This now fails the requested populated-row rule. Move such
  notes outside the sample table or clear that row; do not invent a sample ID.
- Configmaker's descriptor conversion used to trim IDs and lose multi-project
  provenance. The coordinated change preserves the identifiers validated by the
  shared parser, including leading zeros and exact whitespace.
- Configmaker's batch-specific FASTQ discovery had latent argument/name bugs;
  the supported multi-flowcell and `--keep-batch` paths are covered explicitly.
- Repeated samples across submission forms remain supported. Later-form metadata
  updates are reported as warnings, while ambiguous duplicates within one
  worksheet fail. SampleSheet project assignments retain precedence with warnings.

## Deployment order and limits of local checks

Merge/deploy gcf-tools first, then BFQ. Promote the matching dependency to the
branch used by production builds before building the production BFQ image.
The BFQ CI job pins the reviewed companion commit; the test build command above
deliberately uses the development branch so both PRs can be tested before merging.

Automated tests cover parsing, standalone initialization, scheduling order,
curation/refresh, cleanup guards, read-only commands, reporting and notification
identity with fake SMTP. Real demultiplexing, the complete Snakemake analysis,
site mounts/permissions, Docker builds and actual mail delivery require the
server checks above.
