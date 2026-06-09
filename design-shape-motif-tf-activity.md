# Design notes: inferring TF-activity changes from TSS shape + motifs

Dev note (2026-06-08), not user-facing. Follow-up to the question of whether
motif position plus differential TSS shape can infer changes in TF activity at
promoters. Conclusion: yes, but only as a population-level associational
inference. The shape change is the phenotype, the motif assigns a candidate TF
identity, and orthogonal occupancy (ATAC footprint / ChIP) or a perturbation is
what promotes a nomination to a claim. Below are the two analyses, plus the
depth-matching step that must run before either.

## Shared prerequisites

- Per-condition TSS set + quantification, restricted to a testable set with
  `filterByExpression()`.
- Per-condition `tssShape()` output, which already gives `maxPos`, `shannon`,
  `simpson`, and the `pctN` interquantile-width positions per promoter per
  sample.
- New primitive needed: motif scanning. Scan a window around each promoter
  (defined relative to `maxPos`) against a PWM library, recording for each motif
  its presence / count / best-score and its offset and strand relative to the
  dominant TSS. Needs genome sequence (a `BSgenome` or FASTA) plus PWMs
  (JASPAR / CIS-BP via TFBSTools or universalmotif) and `Biostrings::matchPWM`.
  Sketch: `scanPromoterMotifs(TSS, genome, pwms, window = ...)`.

## 0. Depth-matching (runs first, gates both analyses)

- Why: Shannon entropy and IQ-width are downward-biased / distorted at low tag
  counts, so unequal coverage between conditions fakes a shape change.
- Must operate on raw 5' counts, not RPM tracks. `tssShapeFromBigWig` reads
  `.rpm.*.bw`; subsampling needs integer counts, so this has to point at
  un-normalized count tracks instead.
- Method: for each promoter, downsample per-position counts in every compared
  sample to a common target (the minimum count across samples, or a fixed cap),
  recompute shape on the subsampled counts, and average over a few draws to cut
  variance. Drop promoters below a minimum count in any sample, and log how many
  are dropped (no silent truncation).
- Output: depth-matched shape matrices that are a drop-in replacement for the
  raw `tssShape()` output feeding the two analyses.
- Sketch: `tssShapeMatched(TSS, countSignalList, target = NULL, draws = ..., minCount = ...)`.

## 1. MARA-style delta-shape ~ motif regression (population-level)

- Response: a delta-shape vector per promoter from the depth-matched shape, e.g.
  delta-entropy, delta-width, delta-mode-shift, or delta-directionality (one per
  model).
- Predictors: the per-promoter motif matrix (presence / count / best-score),
  optionally position-binned (core window vs upstream) so position enters the
  model.
- Covariates: baseline shape, GC content, coverage, promoter class.
- Model: regularized multivariate regression (glmnet, ridge or elastic-net) to
  deconvolve collinear motifs. The signed coefficient per motif is its
  shape-activity change. Significance by permutation or stability selection,
  since glmnet p-values are awkward.
- Output: ranked motif table with signed coefficients.
- Sketch: `shapeMotifActivity(deltaShape, motifMatrix, covariates, alpha = ...)`.

## 2. Motif-anchored delta-signal meta-profile (position-resolved)

- For a chosen motif, gather all occurrences near promoters with genomic
  position + strand.
- Extract per-position 5' signal in each condition over a fixed window, centered
  and strand-oriented on the motif (flip minus-strand windows; decide
  sense/antisense handling). `genomation::ScoreMatrix` (already a dependency)
  does the windowed extraction.
- delta-signal per position = condition B minus condition A, aggregated across
  occurrences (mean or median + CI), giving a meta-profile of delta-signal vs
  distance-from-motif.
- A localized peak/dip at a fixed offset is the signature of a direct local
  effect. Compare against motif-free or position-shuffled control windows for
  significance.
- Output: a data.frame of distance vs aggregated delta-signal (+ control),
  plottable.
- Sketch: `motifAnchoredDelta(countSignalListA, countSignalListB, motifGR, window = ...)`.

## Packaging

- New file(s): `R/motif.R` (scanning) and `R/shapeActivity.R` (the two analyses
  + depth-matching), or fold into `shape.R`.
- New deps: Biostrings (matchPWM), TFBSTools or universalmotif (PWMs), a
  `BSgenome` or FASTA reader, glmnet (regression). Heavy ones can sit in
  Suggests with runtime `requireNamespace()` checks. Add package-wide
  `@importFrom` in `tssUtils.R` per repo convention; keep camelCase, 2-space
  indent, and the `@author` tag on exports.

## Caveats / validation

- Motif presence is not occupancy. The regularization handles motif collinearity
  but not the presence-vs-binding gap.
- Shape-based inference only sees positioning / architectural TFs. Recruitment-
  only factors show up in level, not shape, so expect fewer and subtler shape
  hits than level hits.
- Core-promoter motifs (TATA at about -28, Inr, DPE, TCT) at stereotyped offsets
  are the most reliable; upstream sequence-specific TFs need the pooled model.
- Confirm any nomination with an ATAC footprint or ChIP in both conditions, or a
  TF knockout / motif mutation (the motif-anchored delta should vanish if the
  link is causal).
- Depth-matching is mandatory in front of both, because a coverage-driven
  entropy or width shift otherwise masquerades as TF activity.
