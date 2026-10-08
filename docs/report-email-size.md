# Analysis report email size

Issue: [#146](https://github.com/gcfntnu/gcf-bfq/issues/146).
BFQ keeps analysis summaries deliverable when report attachments exceed the
message budget. This changes email attachment selection only: reports, FASTQs,
archives, state schema and processing/restart behavior remain unchanged.

## Limit and source

Set the operational budget in `/config/bcl2fastq.ini`:

```ini
[Email]
max_message_bytes = 20000000
```

The value is a positive integer in **bytes**, not raw HTML bytes or MiB.
The default is 20,000,000 bytes (20 MB decimal, approximately 19.07 MiB).
This is a BFQ fallback policy, **not a verified production SMTP limit**.
No production rejection transcript, relay advertisement or separate
per-attachment restriction was supplied for this implementation. Confirm those
on the integration server before merging operational behavior; set the budget
to the verified effective restriction, allowing for downstream relays.

For every analysis-notification attempt BFQ performs EHLO/HELO negotiation before
DATA. A positive numeric `SIZE` advertisement lowers the effective budget when
it is smaller than the configured value. An absent, empty or zero advertisement
does not supply a fixed maximum; an invalid advertisement is ignored with a
warning. The configured value always remains an upper bound. There is no inferred
per-attachment limit and no automatic resend after SMTP rejection.

[RFC 1870](https://www.rfc-editor.org/rfc/rfc1870), sections 3 and 5, defines the
advertised maximum and message size in octets, including headers, bodies and
CRLF, excluding SMTP dot stuffing and the DATA terminator.
[RFC 2045](https://www.rfc-editor.org/rfc/rfc2045), section 6.8, defines Base64
encoding and line wrapping. Raw file length therefore cannot establish that an
email fits. BFQ serializes the complete MIME message using SMTP line endings,
including its actual encoding, omission notices, filenames and envelope-dependent
SMTPUTF8 policy, and measures those bytes.

## Selection and retry policy

If all reports fit, they are attached as before. Otherwise BFQ removes additional
10x/Parse HTML reports first, legacy sequencing HTML attachments next, then
MultiQC reports. Within each priority it removes larger encoded parts first;
ties use the actual report path in ascending order. This order applies across
every project on the flowcell and is deterministic. BFQ stops as soon as the
message fits; it does not re-add lower-priority reports after discarding them.

Every candidate is measured again with all its omission notices. Plain-text and
HTML bodies identify every omitted report, the effective limit and the exact
saved flowcell output path/filename; HTML paths are escaped. A report that cannot
fit with the original summary even as the sole attachment is distinguished from
one omitted due to the combined message budget (including other omission notices).
The summary does not describe omitted MultiQC reports as attached.

All attachments may be omitted. Delivery then counts as successful if SMTP accepts
the summary. If the summary plus all location notices is itself oversized, BFQ
records an actionable notification failure before DATA instead of sending a known
oversized message or silently truncating the summary. Missing/unreadable required
reports remain composition failures; they are not classified as size omissions.

Saved version-1 and version-2 analysis payloads use the same budget policy,
including legacy sequencing attachments. Retries use the current configuration,
recipients and relay advertisement with the saved run context and stable
notification identity. Genuine SMTP rejections remain failures; interruptions
after sending begins remain uncertain. Notification retries do not rerun
processing or modify outputs. Static email configuration changes require the
usual daemon restart; a fresh manager command loads the current configuration.
Early sequencing and finalization/error notifications are outside this policy.

Logs record the configured/advertised/effective byte limits, final serialized
size, each omitted report path, message sizes before/after removal, its size as
the sole report with the summary, and the omission reason. Explicit SMTP response
errors log their code and response; the established notification state retains
delivery failure evidence.

## Verification and remaining server checks

Local regression tests cover encoded-size and exact-byte boundaries, notice
overhead requiring another omission, supplementary/MultiQC priority, multiple
projects and deterministic ties, escaped output paths, summary-only delivery,
invalid settings, relay advertisements, pre-DATA failures, and SMTP rejection.
Mocked daemon/manager tests cover both payload versions, preserved processing and
output bytes/timestamps, an explicit retry after correcting the budget, successful
fallback delivery and suppression of duplicate sends. Local tests send no mail.

For manual integration, use disposable outputs and approved recipients:

1. Record the observed oversize SMTP rejection code/text and the actual relay's
   EHLO `SIZE` advertisement. Ask the relay operator about downstream message
   limits and any per-attachment restriction; do not assume `SIZE` covers those.
2. Configure the verified budget, restart the test instance, and inspect a small
   multi-project analysis notification. All reports should arrive unchanged.
3. Exercise combined oversize with a 10x/Parse flowcell. Verify supplementary
   reports are omitted ahead of MultiQC; check notices and exact existing paths
   in both plain-text and HTML clients, plus logged wire size and effective limit.
4. Exercise an individually oversized MultiQC report and a summary-only message.
   Verify receipt, report preservation, notification `sent` and completed stages.
5. For an existing failed current or legacy analysis notification, correct the
   configured budget and run `fm retry-notifications RUN_ID --kind processed`.
   Verify delivery without processing/output changes. Repeat to verify no duplicate.

The requested integration image is built from the issue checkout with
`bash build-tag-push.sh test dev-test`. The wrapper builds and pushes
`gcfntnu/bfq:dev-test`; it does not start the daemon or deploy the image. Its tools
and workflows use the wrapper's existing `bfq-dev` defaults, which differ from
the immutable local-test companion baseline. Actual mail receipt, image runtime,
scientific outputs and mounted paths remain manual integration evidence.
