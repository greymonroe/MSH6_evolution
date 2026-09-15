# Revision analyses (2026)

Code and result tables for the analyses added during revision of the MSH6 histone-reader
manuscript (eLife Reviewed Preprint RP-RA-2024-105016). The 2024 analyses are unchanged
and live in `code/` and `data/` at the repository root.

Two independent analyses were added:

1. **Exon-cassette architecture** — where the intron falls *inside* the PWWP and Tudor
   reader domains, across every reader-containing protein with a resolvable exon
   structure, used to ask whether the metazoan PWWP was acquired as a two-exon cassette.
2. **Protostome co-evolution** — whether presence/absence of orthologs across 190
   Protostome genomes correlates with presence/absence of the MSH6-PWWP fusion.

---

## 1. Exon cassette — `exon_cassette/`

| script | what it does |
|---|---|
| `fetch_msh6_exon_structure.py` | MSH6 exon structures from NCBI |
| `fetch_exon_data.py`, `fetch_all_exons.py`, `fetch_plant_exons.py` | exon structures for the wider species panels |
| `analyze_exon_cassette.py` | intron positions relative to domain boundaries |
| `pwwp_donor_screen.py` | InterPro PF00855 (PWWP) → UniProt → EMBL → GenPept `coded_by`; computes the PWWP-internal exon break for every PWWP protein |
| `tudor_donor_screen.py` | the same for Tudor domains |
| `ancestral_state.R` | ML marginal ancestral states (`ape::ace`, `type="discrete"`, `model="ARD"`) |
| `make_figure.py`, `make_figure_v2.py`, `build_100sp_figure.py` | exon/domain architecture diagrams |
| `fig_exon_cassette*.R` | manuscript figure panels |
| `PWWP_DONOR_BLAST_2026-09-14.md` | **read this before using the donor-screen numbers** — records the exact BLAST commands, databases, and two corrections made to the Methods |

Key result tables are in `data/revision_2026/exon_cassette/`. The central one is
`pwwp_all_proteins_breaks.tsv` (19,273 PWWP-containing proteins; 10,337 with a resolvable
PWWP-internal exon break, column `break_from_pwwp_start`).

Reproducing the manuscript's numbers from that table:

```python
import csv, re
from collections import Counter
rows = list(csv.DictReader(open("pwwp_all_proteins_breaks.tsv"), delimiter="\t"))
wb   = [r for r in rows if r["break_from_pwwp_start"] not in ("", "NA")]
len(rows), len(wb)                                        # 19273, 10337
sum(1 for r in wb if r["break_from_pwwp_start"] == "62")   # 117

nsd = [r for r in wb
       if "NSD" in ((r["gene"] or "") + " " + (r["product"] or "")).upper()]
len(nsd)                                                   # 403
Counter(r["break_from_pwwp_start"] for r in nsd).most_common(1)   # ('30', 111)
sum(1 for r in nsd if r["break_from_pwwp_start"] == "62")  # 0
```

The +30 modal NSD boundary and the absence of any NSD protein at +62 are robust to how
the NSD family is matched (token match, substring, or including WHSC aliases all give
modal +30, n = 111, and zero at +62).

### On the donor screen — what it can and cannot show

The BLASTp donor screen ranks candidate sources of the MSH6 PWWP, and its top non-MSH6
hits are NSD-family H3K36 methyltransferases. It does **not** establish that NSD is the
donor, and the manuscript does not claim it is: the similarity is modest (46% identity
over 67 of 104 aligned residues) and the NSD family carries a different PWWP-internal
exon boundary.

Pairwise BLASTp is a low-sensitivity instrument *within* this domain family. The same
MSH6 PWWP query recovers only 5,612 of the 15,307 PWWP domains annotated by the InterPro
PF00855 profile HMM — it misses 63% of real PWWPs. So a BLAST non-hit is not evidence of
non-homology. Of the 117 proteins with a +62 break, 90 are MSH6 orthologs and the
remaining 27 are PWWP domains too divergent for pairwise detection from the MSH6 query;
a donor is neither identified nor excluded among them. See
`PWWP_DONOR_BLAST_2026-09-14.md` for the commands and database details.

BLAST databases are not committed (hundreds of MB). They are rebuildable with `makeblastdb`
from the proteomes listed in `PWWP_DONOR_BLAST_2026-09-14.md`; the query sequences
(`query_msh6_*.fasta`) are included.

---

## 2. Protostome co-evolution — `protostome/`

> **Read `CORRECTION_2026-09-14.md` first.** A bug in the topGO setup invalidated the
> original GO-enrichment numbers, and everything downstream of it was re-run. The scripts
> and tables published here are the corrected versions; the correction note records
> exactly what changed and why.

| file | what it is |
|---|---|
| `CORRECTION_2026-09-14.md` | the bug, its effect, and the before/after numbers |
| `orthodb.R` | the corrected analysis — OrthoDB presence/absence, Fisher screen, topGO enrichment, Pagel correlated-evolution tests |
| `fig_protostome_corrected.R` | produces the submitted Figure S7 |

Briefly, the bug: `annFUN.gene2GO` subsets its annotation by name, and the `gene2GO` list
carried duplicate names because OrthoDB v11 nests the same gene at multiple taxonomic
levels. Name-based subsetting returns only the first match, silently dropping ~19,000 of
26,052 annotated orthologs from the universe. Two smaller issues rode along:
`p.adjust()` was left at its default `"holm"` rather than `"BH"`, and one GO term was
duplicated in the results table.

The consequence is a change in interpretation, not just in numbers. After correction,
DNA repair is no longer the top-ranked enriched term (it falls to 8th of 2,124), so the
enrichment is used only to define a candidate set. The primary result is the
BH-corrected Pagel test: of 24 DNA-repair genes with nominal P < 0.05, **12 pass
Benjamini-Hochberg** across all 179 orthologs tested.

Result tables are in `data/revision_2026/protostome/`:
`all_pagel_results.csv` (179 tests, BH-adjusted), `repair_Pagel_sig.csv`,
`repair_candidates_unique.csv`, `TableS5_corrected.csv` (the 24 genes, input to Figure
S7 panel b), and `Protostome_DNArepair_Pagel_sig_by_Name.csv`.

`orthodb.R` expects the OrthoDB v11 inputs under `files/orthoDB/` in the
[`tol_reader_repair`](https://github.com/greymonroe/tol_reader_repair) repository, where
this analysis was originally run; the file is copied here so that the corrected code is
reachable from the repository the paper cites.

### Deliberately not published here

`verify_protostome.R` — written to re-derive the *original* numbers, so it faithfully
reproduces the bug. It was useful in finding the error and is described in
`CORRECTION_2026-09-14.md`, but publishing it next to the corrected code would be
misleading.

The permutation test in the superseded `ROBUSTNESS.md` was computed from the buggy
pipeline, its generating script was never saved, and it was not re-run. The permutation
claim has been removed from the manuscript; it is not load-bearing under the corrected
framing, where the enrichment only selects candidates.
