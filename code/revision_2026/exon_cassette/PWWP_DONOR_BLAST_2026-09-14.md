# PWWP donor BLAST — verification and correction, 2026-09-14

The Methods claim *"BLASTp of the MSH6 PWWP domain against proteomes from 17
diverse eukaryotes; NSD-family H3K36 methyltransferases were the closest non-MSH6
hit (46% identity to human NSD2, E = 2 × 10⁻¹⁴)"* had no saved output — it traced
only to the prose in `R1-7_text_additions.md`. Re-run from the surviving BLAST
databases to verify it.

## Commands

```bash
# 16-proteome screen (the purpose-built database; REPORTED IN THE MANUSCRIPT)
blastp -query query_msh6_pwwp.fasta -db combined_20sp_db \
       -outfmt "6 sseqid pident length evalue bitscore stitle" \
       -evalue 10 -max_target_seqs 5000 -num_threads 6 \
       -out pwwp_donor_blast_16proteomes.tsv

# 6-proteome screen (the smaller database; kept for provenance)
blastp -query query_msh6_pwwp.fasta -db combined_proteomes_db \
       -outfmt "6 sseqid pident length evalue bitscore stitle" \
       -evalue 10 -max_target_seqs 5000 -num_threads 6 \
       -out pwwp_donor_blast_6proteomes.tsv
```

Query = `query_msh6_pwwp.fasta`, human MSH6 PWWP (NP_000170.1), 104 aa.

## Result — the biological claim holds

The closest non-MSH6 hit is NSD-family H3K36 methyltransferase in both databases,
followed by NSD3, then zebrafish NSD2. **Percent identity is 46.269% either way**
(identity does not depend on database size). The conclusion — NSD family is the
closest extant relative — is verified.

## Two errors in the Methods as written, both now corrected

**1. Species count.** `combined_20sp_db` contains **16 distinct species**, not 17
(Amphimedon queenslandica, Arabidopsis thaliana, Caenorhabditis elegans, Danio
rerio, Dictyostelium discoideum, Drosophila melanogaster, Gallus gallus, Homo
sapiens, Nematostella vectensis, Neurospora crassa, Oryza sativa, Plasmodium
falciparum, Saccharomyces cerevisiae, Schizosaccharomyces pombe, Tetrahymena
thermophila, Thalassiosira pseudonana). The "17" came from counting organism
*labels*, where Dictyostelium appears twice ("Dictyostelium discoideum" and
"Dictyostelium discoideum AX4" — same species). The db name "20sp" reflects an
intended, not achieved, species count.

**2. E-value belonged to the wrong database.** BLAST E-values scale with database
size, so the two runs do not agree:

| database | residues | best human NSD2 hit | best non-MSH6 hit overall |
|---|---|---|---|
| `combined_proteomes_db` (6 species) | 191,701,571 | **2.40e-14** | Homo NSD2, 2.40e-14 |
| `combined_20sp_db` (16 species) | 335,139,280 | **4.20e-14** | Gallus NSD2, 2.57e-14 |

The manuscript reported "17 diverse eukaryotes" (implying the larger screen) but
quoted **2 × 10⁻¹⁴**, which is the *6-proteome* value. Corrected to report the
16-proteome screen throughout: **46% identity to human NSD2, E = 4 × 10⁻¹⁴**.

## Related exon-boundary numbers, independently re-verified against `pwwp_all_proteins_breaks.tsv`

- 10,337 PWWP-containing proteins have a resolvable PWWP-internal exon break
  (of 19,273 rows total). **Matches the manuscript.**
- NSD family: 403 proteins with a break, modal boundary **+30 aa (n = 111)**.
  **Matches the manuscript exactly.** (Widening the match to include WHSC aliases
  gives 416 with a break; the mode is +30, n = 111, either way.)
- **Zero NSD-family proteins have a break at +62**, supporting "NSD ... do not
  share the MSH6-characteristic exon boundary."
- A +62 break occurs in **117** proteins. BLASTing the human MSH6 PWWP against
  `pwwp_db` (E < 1e-3) shows **90 of 117 are MSH6 orthologs**; the remaining 27
  have no significant similarity to the MSH6 PWWP domain at all (mostly fungi and
  insects). So no protein combines close MSH6-PWWP sequence similarity *with* the
  +62 architecture — which is what the corrected Methods sentence now says. The
  earlier, stronger claim ("no extant PWWP-containing protein shares the +62
  boundary") was false as literally written and had already been replaced.
