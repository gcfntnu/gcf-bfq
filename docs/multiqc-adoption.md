# MultiQC adoption and evaluation

This is the BFQ handoff for [issue #173](https://github.com/gcfntnu/gcf-bfq/issues/173)
and the companion [gcf-workflows PR #222](https://github.com/gcfntnu/gcf-workflows/pull/222)
([issue #221](https://github.com/gcfntnu/gcf-workflows/issues/221)).
The candidate replaces the old MultiQC fork with upstream MultiQC and the
[GCF extension](https://github.com/gcfntnu/gcf-multiqc/pull/1).
Image packaging and the coordinated source revisions are tracked in
[gcf-docker #1](https://github.com/gcfntnu/gcf-docker/pull/1).
This checkpoint prepares reviewable source changes; it does not approve deployment.

The published candidate tag is `gcfntnu/gcf-multiqc:1.35-gcf0.1.0-rc2`. Its immutable
reference, selected by the companion workflow change, is:

```text
gcfntnu/gcf-multiqc@sha256:8c53227c570280641c51738a8e7308ae11cddd00a1e0cd1092e5cef7d92cc9e9
```

## Configuration and compatibility

MultiQC image selection remains workflow-owned, with defaults in
`gcf-workflows/docker.config`. BFQ has no MultiQC version setting in its INI
files, and its Python environment must not install the reporting engine: that
would replace the BFQ `multiqc` wrapper.

| Caller | Image selection | Effect of adopting the new workflow checkout |
| --- | --- | --- |
| BFQ installed `multiqc` wrapper, including sequencing QC | `/opt/gcf-workflows/docker.config`, or the explicit `GCF_WORKFLOWS_DOCKER_CONFIG` override | The next invocation uses the selected map. The wrapper prints the resolved image on stderr. |
| Project and aggregate analysis reports run by Snakemake | Explicit project `docker.multiqc` configuration wins; missing keys are filled from that analysis workdir's `src/gcf-workflows/docker.config` | New workdirs receive the workflow revision shipped with BFQ. Retained workdirs keep their own copied mapping and workflow code; explicit project overrides still win. |

Changing the wrapper override does not update retained workflow copies or project
overrides. Record the actual mapping and workflow revision used for each
evaluation. Keep any review changes in a copied evaluation workdir; do not
rewrite operational state, delete outputs or alter FASTQs to force report
regeneration.

The sequencing-report command passes the BFQ-owned YAML configuration,
`--force`, `--no-megaqc-upload`, output directory and filename, with explicit
`custom_content`, `interop` and the detected `bclconvert` or `bcl2fastq` module.
Its HTML summary is supplied as
an input file; the sample assignment table is in the owned configuration.
Separate conversion passes select only custom content and InterOp. Previously
this path repeated `custom_content` in the module order, causing duplicate BFQ
sections under MultiQC 1.35 and failure in strict mode. The adoption change
selects each module once, using the same list for the configuration and CLI.
Native bcl-convert/bcl2fastq selection and combined-pass count semantics remain
unchanged. Report paths, failure handling, zero/unavailable sample presentation and
sequencing-notification scheduling remain as documented in the
[operator guide](../README.md#early-sequencing-qc-and-the-three-notifications).

## Local verification at this checkpoint

Run the normal checkout-local development gates:

```console
python3.11 scripts/dev.py setup
python3.11 scripts/dev.py check fast
python3.11 scripts/dev.py check all
```

`scripts/check_install.py` now invokes the actual installed `multiqc` command
with a fake Apptainer executable. The check verifies the candidate digest,
mapping replacement, paths and arguments containing spaces, inherited stdin,
stdout/stderr, nonzero exit propagation and failure before execution when the
mapping is missing. The fake child must inherit the SMTP guard. `fast` checks
the editable install; `all` repeats this check against the clean installed
wheel. Sequencing-QC regression tests cover the separate-pass selector fix,
native selectors, preserved zero/unavailable counts and existing failure paths.

These checks do not install MultiQC, execute containers, render scientific
reports, send mail or run the facility integration setup. Passing them establishes
the BFQ command boundary, not the quality of reports rendered by the candidate.

## Later build and evaluation

Build only after the companion PRs have reached the agreed review checkpoint.
The BFQ test image must include the reviewed workflow revision; building BFQ
against unchanged `bfq-dev` workflows would keep the previous image mapping.
For example, from the clean reviewed BFQ branch, set
`GCF_WORKFLOWS_REVIEW_REV` to the full reviewed commit from
`issues/221-multiqc-adoption`, then build a local image:

```bash
: "${GCF_WORKFLOWS_REVIEW_REV:?Set the full reviewed workflows commit}"
[[ "$GCF_WORKFLOWS_REVIEW_REV" =~ ^[0-9a-f]{40}$ ]] || exit 1
test -z "$(git status --porcelain)" || exit 1
BFQ_REVIEW_REV=$(git rev-parse HEAD)
docker build -f dockerfile-test -t bfq:multiqc-adoption-review \
  --build-arg BASE_IMAGE=gcfntnu/bfq:base-261012 \
  --build-arg BFQ_REVISION="$BFQ_REVIEW_REV" \
  --build-arg GCF_WORKFLOWS_BRANCH=issues/221-multiqc-adoption \
  --build-arg GCF_WORKFLOWS_REV="$GCF_WORKFLOWS_REVIEW_REV" \
  --build-arg GCF_TOOLS_BRANCH=bfq-dev \
  --build-arg GCF_TOOLS_REV=1daa9822084bc3a3df4b94102c13c6352fa57ecd .
```

The tools revision above is the BFQ development dependency baseline, with no
local override. Record the BFQ/workflow/tools commits and built image labels in
the evaluation record. `build-tag-push.sh test` is a separate **build and publish**
action: it pushes immediately and resolves the selected companion branches when
invoked. Its legacy `-q` MultiQC branch option is unused. Do not use that helper
for a build-only review, or omit `-w issues/221-multiqc-adoption` if using it for
an explicitly agreed later publication of this candidate.

Evaluate reports from copied input logs and configurations in a separate local
area. Keep the old reports and their data exports as references; keep copied
operational inputs out of Git, published images and public review attachments.
Do not run a processing daemon or trigger real notifications for report comparison.

1. Verify the built image contains the reviewed workflow map and the BFQ wrapper.
   Confirm Apptainer can resolve the pinned image under the intended mounts and
   user permissions, and capture its version and selected image reference.
2. Render sequencing QC from copied bcl-convert and bcl2fastq inputs. Check the
   BFQ summary, expected/zero/unavailable sample assignment table, native
   demultiplexing sections, optional InterOp behavior, HTML/data outputs and
   failure reporting. Compare sample coverage and units, not byte-identical HTML.
3. Exercise the companion workflow's project and aggregate report invocations
   on copied report inputs. Include Parse reports and the GCF custom content;
   QIIME2 and old mkfastq are deprecated and STARsolo is low priority. Review
   documented metric corrections in the extension before treating differences
   from the fork as regressions.
4. Check a copied retained workdir separately from a fresh one. Updating the
   BFQ image does not rewrite its `src/gcf-workflows` snapshot. A reporting-only
   BFQ rerun preserves project analysis reports, so it is not evidence that
   those reports were regenerated with the candidate. Use the explicit report
   invocation in the copied workdir when comparing existing analysis logs.

Operational report comparison, container/mount behavior and facility integration
remain manual acceptance gates before merge and deployment. The BFQ development
checks do not execute those acceptance tests.

## Rollback

Keep the previous BFQ image identity, workflow revision and image mapping with
the evaluation record. Re-selecting the previous workflow map changes the BFQ
wrapper's next invocation; an explicit wrapper override must also be restored
or removed deliberately. Analysis reports use their own retained workflow
copies, so rollback must account for the exact copied mapping used there as
well. Do not infer their selection from the image shipped at `/opt/gcf-workflows`.
Explicit project `docker.multiqc` overrides must also be restored when applicable.

Preserve the successful old reports and evaluation outputs for comparison.
Rollback does not require changing state JSON, deleting successful reports,
rerunning demultiplexing, or modifying FASTQs and their valid MD5 manifests.
Any later operational retry follows the existing BFQ lifecycle and lease rules.
