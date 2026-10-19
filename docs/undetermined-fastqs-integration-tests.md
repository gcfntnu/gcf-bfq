# Undetermined retention: manual integration checks

Use disposable test conversions and a test mail environment under maintainer
control. These checks are required before merging issue #157; the local suite
uses synthetic trees and doubles for converters, scientific tools and SMTP.
Keep actual facility identifiers, storage details and raw logs private; publish
only a synthetic/general result summary.

## Naming and ownership

The policy recognizes only the converters' reserved sample-zero filenames:
`Undetermined_S0[_Lnnn]_<R1|R2|R3|I1|I2>_nnn.fastq.gz`. Lane and file-number
fields contain three digits; matching is case-sensitive. Ordinary samples use
positive sample numbers, including samples named `Undetermined`. BFQ leaves
reserved names intact during renaming. Previously renamed files lacking the
reserved S0 identity are not guessed at or automatically deleted.

The naming contract comes from the
[BCL Convert output guide](https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/OutputFiles_swBCL.htm),
[bcl2fastq 2.20 guide](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf),
[Cell Ranger mkfastq output examples](https://www.10xgenomics.com/support/software/cell-ranger/7.2/analysis/inputs/cr-mkfastq)
and [Cell Ranger ATAC mkfastq output examples](https://www.10xgenomics.com/support/software/cell-ranger-atac/2.1/tutorials/inputs/generating-fastqs-mkfastq).
The 7-Zip recursive masks retain the reserved name/read identity and field
widths; its single-character wildcard cannot enforce numeric fields. They also
filter files reached through archive directory aliases, without scanning those
external trees during cleanup.

## Before merge

- Exercise BCL Convert without lane splitting, bcl2fastq with lane splitting,
  and the supported Cell Ranger, Space Ranger and Cell Ranger ATAC mkfastq
  commands/versions available in the integration environment. Inspect actual
  output and staging layouts (including `outs/fastq_path`, Martian directories,
  aliases, index reads and ATAC R3). Confirm every generated unassigned FASTQ
  retains the reserved identity and is removed before scientific analysis.
- Compare conversion statistics before/after cleanup: assigned and Undetermined
  counts/percentages, unknown-index summaries, InterOp, native reports and logs
  remain available. Regenerate early QC from retained statistics and inspect
  actual MultiQC output. Sample FASTQs and valid MD5 manifests must remain
  unchanged, including similarly named samples and multiple project directories.
- List new delivery archives with the deployed 7-Zip version. Check both normal
  project inputs and the optional nested mkfastq directory/alias: no recognized
  Undetermined read/index FASTQ should appear; metadata, statistics, reports,
  project reads and checksums remain included. Existing completed archives and
  their checksums must remain unchanged until an explicit finalization rerun.
- As the unprivileged test service account, inject an unlink permission failure
  after one removal, then separately interrupt cleanup. Verify completed
  demultiplexing, an analysis failure/interruption and useful path/count/byte
  logs. Fix the cause and use `fm rerun RUN_ID --from analysis`; confirm no
  converter rerun or duplicate sequencing email, and successful retry with
  already-absent files. Repeat cleanup/recovery without restoring these reads.
- Include file/dangling symlinks and external directory aliases in the disposable
  tree. Only selected file links are unlinked; external targets and directories
  remain intact. A failed conversion must retain diagnostic FASTQs. A test SMTP
  failure or failed early report must not prevent cleanup or completed processing;
  use the normal notification/report retry procedures afterward.

No scientific QC acceptance threshold, completion marker, configuration switch,
companion-repository change or multi-project orchestration is introduced.
