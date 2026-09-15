# Correction, 2026-09-14 — read this before trusting anything else in this directory

A bug in the topGO setup in `code/orthodb.R` (repo `tol_reader_repair`) invalidated the
GO enrichment numbers produced by the original run. Everything downstream of that
enrichment was re-run. Files in this directory are now split into two groups.

## The bug

`annFUN.gene2GO` subsets its annotation with
`gene2GO[intersect(names(gene2GO), feasibleGenes)]`. The `gene2GO` list carried
**duplicate names** — OrthoDB v11 nests the same gene at multiple taxonomic levels
(at2759 Eukaryota, at33208 Metazoa, ...), up to 12 copies of one gene. List subsetting by
name returns only the **first** match, so ~19,000 of 26,052 annotated orthologs were
silently dropped from the universe and the enrichment was computed on an arbitrary
one-per-name subset.

Two smaller issues rode along: `p.adjust()` was left at its default `"holm"` rather than
`"BH"`, and `GO:0006334` appeared twice in the results table.

## What changed

| | before (buggy) | after (corrected) |
|---|---|---|
| topGO universe | ~7,278 orthologs | 26,052 annotated / 24,766 assigned to a testable term |
| Fisher-screened ortholog groups | reported as "26,052" | **67,241** (all groups present in >= 20 species) |
| GO:0006281 rank (classic Fisher) | 1 | **8** of 2,124 terms |
| GO:0006281 p | 9.9e-20 | **4.5e-11** (90 of 883 annotated significant vs 43.4 expected) |
| top term | DNA repair | protein ubiquitination (2.5e-13) |
| Pagel multiple-testing correction | Holm | **Benjamini-Hochberg**, across all 179 tests |
| DNA-repair genes co-evolving with MSH6-PWWP | "24 significant" | 24 nominal, **12 pass BH** |

Because the enrichment no longer ranks first, the manuscript no longer treats it as
evidence of co-evolution. It now defines a candidate set only; the BH-corrected Pagel
correlated-evolution test is the primary result.

## Current files (regenerated 2026-09-14 — these are authoritative)

- `all_pagel_results.csv` — all 179 Pagel tests (90 GO:0006281 + 89 comparison-set), BH-adjusted
- `repair_Pagel_sig.csv` — the 47 repair ortholog rows with nominal Pagel P < 0.05
- `repair_candidates_unique.csv` — the 90 GO:0006281 candidate orthologs. Column names were
  made explicit; the old `p` / `p.adjust` pair silently mixed Fisher and Pagel statistics.
- `Protostome_DNArepair_Pagel_sig_by_Name.csv` — 24 genes, collapsed across taxonomic
  levels. 1:1 with `revision/TableS5.xlsx` (verified field-by-field).
- `FigS_protostome_coevolution.pdf` / `.png` — regenerated; panel a is classic Fisher,
  panel b colours the 12 BH-passing genes green.

## SUPERSEDED — do not cite, do not copy numbers out of

- `ROBUSTNESS.md` §2b — the permutation test. Its observed statistic (7.2e-19, described as
  matching the saved 9.9e-20) came from the buggy pipeline, so its null is not comparable to
  the corrected analysis. **The permutation claim has been removed from the manuscript and
  the response to reviewers.** It was not re-run: under the corrected framing the enrichment
  is only used for candidate selection, so the permutation is no longer load-bearing.
  Separately, the script that generated it was never saved — it is not reproducible as-is.
- `permutation_null_distribution.csv`, `permutation_null_distribution_1000.csv`,
  `perm_topgo_results.csv` — outputs of that same permutation.
- `Protostome_DNArepair_Fisher_sig_by_Name.csv` — derived from the buggy candidate set.
- `fig_protostome.R` — the pre-correction figure script. It reads the buggy
  `GO_enrich_protostomes.csv`, has no BH colouring, and points at a `hub` path that no
  longer exists. **Use `fig_protostome_corrected.R`**, which is the script that actually
  produced the submitted `FigS7_graphic.pdf` (rescued 2026-09-14 from `/tmp/figS7_new.R`,
  where it had been written and left; its `/tmp/TableS5_corrected.csv` input was rescued to
  `TableS5_corrected.csv` here and the script repointed at it).
- `verify_protostome.R` — written to re-derive the *original* numbers, so it **reproduces
  the bug** (`annFUN.gene2GO` on duplicate-named `gene2GO`, `p.adjust` left at `holm`).
  Kept only as the record of how the bug was found. Not published to the code repository.
- `EXPLAINER.md`, `HOW_IT_WORKS.md`, `MANUSCRIPT_ADDITIONS.md`, `to_apply_protostome.md` —
  written against the pre-correction numbers. Kept as a record of what was originally
  drafted; the manuscript and response to reviewers are the current text.

Corrected source: `code/orthodb.R` in `tol_reader_repair` (5 edits) and
`revision/figures/make_FigS7.tex` for the standalone portal export.
