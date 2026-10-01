# Recording the BFQ 2 baseline after production promotion

Use one permanent annotated Git tag, **`bfq2-baseline`**, to identify the
production promotion commit that establishes BFQ 2. The accompanying
[operational milestone](bfq2.md) documents the changes; the
[tag message](bfq2-tag-message.txt) provides a short explanation in Git itself.

This is a one-time historical marker. It introduces no release schedule,
bundled-release requirement, or automatic major/minor/patch policy. Routine
updates continue through `bfq-dev` to `master`, while BFQ remains version `2`
and production image builds increment `prod2-N` independently.

## Prepare and promote

1. Merge the milestone documentation PR into `bfq-dev` before the production
   promotion, so the tagged tree contains the document and README link.
2. Complete the intended integration checks and review the normal promotion PR
   from `bfq-dev` to `master`. Include the operational milestone link in its
   description, alongside the testing actually performed.
3. Merge that PR with a regular merge commit. Record its full commit SHA from
   `master`; this is the tag target. Do not tag the documentation commit or the
   pre-promotion `bfq-dev` tip as a substitute.
4. When the branches are otherwise unchanged, fast-forward `bfq-dev` to the
   promotion merge commit so both long-lived branches start the next cycle from
   the same history. If new development has already landed, review that history
   rather than resetting or force-pushing either branch. Branch alignment does
   not change which production commit receives the milestone tag.

No promotion SHA or production image number is filled in by these documents:
those facts must come from the actual promotion and build. The Git tag records
the source boundary; the production Docker tag identifies a separately built
environment. Tag creation neither builds nor deploys that environment.

## Create and publish the annotated tag

From a checkout containing these documents, fetch the current remote refs and
replace the placeholder with the full, verified promotion merge SHA:

```bash
git fetch origin master bfq-dev --tags
PROMOTION_COMMIT=REPLACE_WITH_FULL_PROMOTION_MERGE_SHA
git show --no-patch --format=fuller "$PROMOTION_COMMIT"
git merge-base --is-ancestor "$PROMOTION_COMMIT" origin/master
git show "$PROMOTION_COMMIT:docs/bfq2.md"
git tag --list bfq2-baseline
```

Confirm that this is the merged production PR, that the ancestry check succeeds,
and that the milestone document is present at that commit. Review the linked tag
message before creating it. If `bfq2-baseline` already exists, inspect its
annotation and target and stop; do not move, delete, or recreate a published
baseline tag. An existing correct tag makes this step unnecessary.

```bash
git tag -a bfq2-baseline "$PROMOTION_COMMIT" -F docs/bfq2-tag-message.txt
git show --no-patch bfq2-baseline
git rev-parse 'bfq2-baseline^{commit}'
```

Verify that the resolved commit equals the recorded promotion SHA. Then publish
only this tag and inspect the remote result:

```bash
git push origin refs/tags/bfq2-baseline
git ls-remote --tags origin 'refs/tags/bfq2-baseline*'
```

The annotated tag object has its own SHA; the peeled `refs/tags/bfq2-baseline^{}`
entry identifies the commit and should equal the promotion SHA.

## What visitors and automation should use

- **Current development:** clone the default `bfq-dev` branch.
- **Current stable production source:** select `master` explicitly.
- **The historical BFQ 2 starting point:** select `bfq2-baseline` explicitly.
- **A particular built production environment:** select its published `prod2-N`
  image tag; the Git milestone does not freeze companion repository branches,
  external images, or future build inputs.

The tag will remain at the baseline even as later fixes reach `master`.
Automation intended to follow current production should therefore select
`master`, not infer it by choosing the newest tag. Existing tag-following
automation should be reviewed before publication; an explanatory annotation
cannot change how a client selects revisions.

A GitHub Release is unnecessary for this plan. If a public milestone page is
wanted later, it can refer to the same tag and operational document and state
clearly that it is a historical baseline, not the current stable checkout.
Publishing such a page should be a separate deliberate action, without
introducing bundled deliverables or a recurring release process.
