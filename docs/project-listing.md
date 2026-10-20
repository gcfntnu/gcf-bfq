# Project rows in list and search

`fm` and `flowcell-manager` are the same installed command. Their `list` and
`search` commands project recorded state into **one project on one flowcell**
per row. This intentionally replaces the former comma-joined flowcell table.
`list-processed` still shows the compatibility inventory unchanged. This projection
is the foundation for #161; it does not introduce holds or change scheduling.

## Columns and identity

| Column | Meaning |
| --- | --- |
| `project`, `run_id` | Row identity; run IDs normalize historical path spellings. Unknown project identities use `<unidentified>`. |
| `status`, `stage` | Current project outcome and active/next boundary, using the rules below. Filters use these columns. |
| `analysis`, `delivery` | Separate current analysis and finalization outcomes. Delivery describes committed processing, independently of notification acceptance. |
| `last_delivery` | Timestamp of recorded successful delivery, even if a new correction is pending/failed. This does not check files or guarantee that archived outputs remain online. |
| `priority`, `priority_source` | Recorded decimal priority and `samplesheet` (explicit) or `default` provenance. Missing records are `unavailable`, never an invented default. |
| `parent_status`, `parent_stage` | Parent flowcell values, labelled separately from project progress. Inventory-only history uses its historical status and stage `legacy`. |
| `origin`, `output_path`, `archived` | State origin, recorded output path and archive timestamp (`0` when absent), retained for historical visibility. |
| `analysis_generation`, `last_delivery_generation` | Current analysis generation and last successful delivery generation. Historical missing identities are `unavailable`. |

Rows sort lexically by run ID, then project. Repeated inventory entries collapse
by project and normalized run ID, using the lexically first recorded path and
latest nonzero timestamps. A project on two runs remains two rows. State owns
an entire run before filters apply; its stale inventory entries cannot reappear.
Project identities come from recorded project/progress/completion/scheduling data,
without opening SampleSheets, workbooks or output directories.

## Current project mapping

Project progress/commits take precedence over parent flowcell outcomes. Read-only
listing reports recorded status; it does not probe/recover abandoned execution
leases. Use the existing recovery commands to persist an interruption.

| Evidence | Project status/stage | Analysis and delivery |
| --- | --- | --- |
| Current project delivery committed | `completed` / `finalization` (`archived` when recorded archived) | Both completed, even if another project is running/failed. |
| Current project analysis running or failed | Its own outcome / `analysis` | Own analysis outcome; delivery pending. An affected project in a recorded interrupted stage displays interrupted. |
| Current project finalization running or failed | Its own outcome / `finalization` | Completed analysis; own delivery outcome. A recorded interrupted stage displays interrupted. |
| Completed analysis awaiting downstream work | `pending` / `reporting`, or `finalization` after reporting completes | Analysis completed; delivery pending. It never inherits a sibling's archive failure/running state. |
| Shared reporting executing/failed/interrupted | Shared outcome / `reporting` for participating selected, undelivered analyzed projects | Analysis completed; delivery pending. Already delivered projects keep their own completion. |
| New preparing/queued project work | `preparing` or `queued` / next boundary | Pending work, with completed earlier boundaries retained. |
| Shared preflight/demultiplexing | Shared parent outcome / `demultiplexing`; analysis preflight failures use `analysis` | Analysis/delivery pending for participating work. Retained unselected project success remains separate. |
| Not-yet-started project in a recorded project-aware attempt | `pending` / `analysis` | Analysis/delivery pending; another project's failure does not assign it a failure. |
| Historical state without individual progress | Completed shared stages support known completion; ambiguous running/failure detail displays `unavailable` at its boundary | Do not assign the parent's analysis/archive failure to every project. |
| Inventory-only history | `completed` or `archived` / `legacy` | Analysis/delivery/generation/priority detail unavailable; a nonzero inventory timestamp records historical completed processing. |

A selected new scientific attempt invalidates its earlier current completion;
its prior delivery stays visible in `last_delivery`. Current analysis and delivery
must refer to the same generation to establish current completion. Unselected
projects keep their own successful completion and producing priority across
scoped corrections. Preparing state can temporarily retain preceding shared-stage
metadata; invalidated stages are not current completion evidence.

Priority comes first from the project's recorded producing analysis, then the
applicable frozen attempt/completion plan. A queued new analysis has no newly
resolved priority until its plan is frozen. Listing never silently resolves it
from today's mutable inputs. Historical missing scheduling data remains
unavailable, including unselected projects without a producing priority.

## Search, filters and read-only behavior

Search matches case-insensitive literal substrings separately against `project`
and `run_id`. A project-only query returns its matching rows; a run query returns
its project rows subject to filters. Searching across formerly joined project
boundaries is impossible. Arbitrary output path components are not searched.

`--status` and `--stage` match exact project row values, independently of parent
status/stage. `--status legacy` is a compatibility selector for inventory-only
rows; their stage is `legacy` and their historical status is completed/archived.
Unknown selectors return no rows. Both filters may be combined.

Exit conventions stay unchanged: 0 for a successful command including no
matches (`No matching flowcells.`), 2 for invalid arguments/blank search queries,
and 1 for state errors. Listing/search do not create manager directories, recover
state, refresh inputs, alter notifications/completions, inspect outputs or send
mail. Missing/archived output directories do not prevent known history appearing.

Run the [development gate](development.md) and
[server listing checks](flowcell-search-integration-tests.md) when changing this
contract. Local regression tests use synthetic metadata and mock SMTP in children.
