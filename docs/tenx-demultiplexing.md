# 10x conversion and legacy reruns

10x demultiplexing is selected independently of analysis. BFQ reads the fixed
workflow-owned `/opt/gcf-workflows/demultiplexing.config` and image mapping in
`docker.config`. `libprep.config` still selects each project's downstream workflow,
chemistry and analysis defaults. Installing BCL Convert does not upgrade analysis.
The [shared compatibility matrix](https://github.com/gcfntnu/gcf-workflows/blob/bfq-dev/docs/tenx-demultiplexing.md)
is authoritative for scoped assays, retained historical routes and integration evidence.

The initial direct profiles cover ordinary 3-prime v3 GEX and a separate Multiome
GEX-only route. ATAC, joint Multiome, historical v1/v2 and Visium do not inherit the
ordinary GEX read geometry. Unlisted assays fail before conversion; install a
reviewed workflow-owned profile rather than changing the kit to bypass validation.

## Explicit selection and indexes

Use `[CustomOptions]` in each effective project SampleSheet:

```csv
[CustomOptions]
Libprep,10X Genomics Chromium Single Cell 3p GEM Library & Gel Bead Kit v3
Demultiplexer,bcl-convert
Index2Orientation,reverse-complement
[Data]
Sample_ID,Sample_Project,index
library-A,GCF-2099-001,SI-TT-A1
```

Choose the documented instrument/reagent orientation. Named TT indexes require
`Index2Orientation,forward` or `reverse-complement`; BFQ does not infer it from an
instrument name. Literal i7/i5 sequences are used as written. A GA name expands to
all four 8-base i7 oligos with the same Sample_ID; one submission-form row remains
sufficient. A dual name expands to its correct paired i7/i5 assignment. Tools owns
the [offline catalogue and validation](https://github.com/gcfntnu/gcf-tools/blob/bfq-dev/docs/tenx-indexes.md).
Unknown names, incomplete/conflicting pairs and collisions fail before conversion.

Select `Demultiplexer,mkfastq` to use the profile's pinned conversion-only image.
BFQ invokes `bfq-container cellranger-mkfastq mkfastq ...`, or the ATAC/Space Ranger
counterpart. Those image keys are separate from downstream analysis images.
Historical profiles default to this route; unsupported direct conversion fails
with a remedy. There is no automatic retry using another converter or orientation.
The existing non-10x path and `FORCE_BCL2FASTQ` behavior remain unchanged; multi-pair
conversion still rejects that environment override in favor of a structured plan.
Single-sheet 10x routing continues to ignore `FORCE_BCL2FASTQ`, as before; select
its backend through `Demultiplexer` and the workflow-owned profile.

RunInfo.xml must define R1/i7/[i5]/R2 for direct 10x conversion. Every cycle is
allocated, excess index cycles are ignored explicitly, and R1 retains required
cell barcode/UMI cycles. Existing adapter trimming is rejected; it is never
silently enabled on barcode reads. Explicit cycle masks, mismatch allowances and
trimming/read-retention settings are checked against the selected profile.
See the matrix for exact minima and BCL Convert 4.2.4 settings.

`fm rerun RUN --toggle-index1` and `--toggle-index2` continue to be deliberate,
source-relative sequence corrections. For symbolic inputs, shared tools first
expands sequences, retaining symbolic labels in annotations and the original
sheet in the existing reversible input history. Complementing actual sequences
then happens once. Subsequent conversion treats those explicit cells as written.
No symbolic name is reverse-complemented. Confirmation/validation and execution
leases still protect the edit; `--force` does not bypass them.

## Projects, outputs and recovery

Compatible projects use one combined pass by default. Differing supported
per-sample settings are represented in the derived v2 SampleSheet. When global
BCL Convert settings cannot share a pass, assign explicit
`DemultiplexingPass` labels in the affected project sheets. The planner fails
before processing when settings in the same pass disagree. Whole-flowcell
collision and identity validation still applies across passes. Mixed legacy/direct
backend plans fail before execution because their shared assay/index collision
contract is unverified; use a supported modern BCL Convert profile or retain the
standalone historical route.

Separate passes retain reports under `Reports/passes/<label>/`. Project FASTQs are
published without overwriting another pass's output. Early QC combines assigned
libraries once, retains each pass's native statistics and leaves combined
Undetermined/unknown counts unavailable rather than adding overlapping counts.
The existing successful-run Undetermined retention policy applies to all layouts.

10x standard FASTQ names are preserved, including sample-number/lane suffixes.
Configmaker discovers every oligo file for the logical sample, and the consumer
reads its raw R1/R2 files without requiring mkfastq-only directories. The optional
legacy flowcell directory remains an optional archive input. Checksums, scoped
corrections, analysis/reporting/finalization restart boundaries, notification
failure handling and successful snapshot retention use the existing lifecycle.

`provenance/demultiplexing/` retains original sheets and forms, derived converter
sheets, physical RunInfo hash, resolved sequences and row mappings, catalogue
identity, profile bytes/hash, effective settings, selected backend, command and
software/image identity. Each project receives its own delivery provenance;
downstream reruns retain its original conversion identity. Source inputs remain
curated and authoritative; derived sheets do not replace them.
Analysis uses the effective biological read lengths from that retained conversion,
including explicit cycle crops; later metadata edits cannot change retained FASTQs.
Preparing and freezing a single-sheet input pair carries RunInfo with it, so clean
initialization and retries validate the same physical geometry after copying inputs.

## Verification

Run the normal `scripts/dev.py setup`, `check fast` and `check all` with the
coordinated tools baseline. The suite covers named and explicit inputs, four-oligo
identity, cycle/orientation/trimming errors, reversible corrections, failed conversion
without fallback, native names, pass isolation, QC and existing recovery/delivery.
Local execution and SMTP boundaries are mocked. The tiny #143 Snakemake text DAG
remains the only local real-workflow exception.

Authorized server validation uses disposable copied fixtures. Compare assignments,
read counts, sequences and qualities with the established route, then inspect the
intended consumer's barcode retention and outputs, sequencing/analysis QC, checksum
manifests, archives, retained snapshots, real test SMTP and reruns. Compression,
headers and ordering may differ; compressed FASTQ MD5 equality is not the equivalence
criterion. Scientific acceptance and any unverified historical cases are recorded
in the shared matrix. A vendor tiny conversion example is not scientific evidence.
