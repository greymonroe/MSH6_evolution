## Corrected Figure S7 (Protostome co-evolution), 2026-09-14.
## Rescued verbatim from /tmp/figS7_new.R -- the script that produced the
## submitted FigS7_graphic.pdf -- with its /tmp/TableS5_corrected.csv input
## repointed to the copy committed alongside this script.
## Supersedes fig_protostome.R, which was written against the buggy topGO run.
suppressMessages({library(data.table); library(topGO); library(ggplot2); library(patchwork)})
setwd("/Users/greymonroe/repos/tol_reader_repair")
out <- "/Users/greymonroe/Dropbox/Research/MSH6_tudor/eLife_resubmission_2026"
theme_set(theme_classic(base_size = 6))
Multiple_threshold_test <- function(allScore) allScore == 1

## ---------- rebuild the corrected topGO result ----------
wn  <- names(fread("files/orthoDB/OrthoDB_species_wide.txt", nrows = 0))
X   <- fread("files/orthoDB/odb11v0_OG_xrefs.tab.gz"); colnames(X) <- c("Ortho_ID","exDB","exDB_ID","N"); X <- X[Ortho_ID %in% wn]
ch  <- fread("files/orthoDB/chitests_protostomes.csv")
bio <- merge(ch, X, by="Ortho_ID")[exDB=="biological_process"]
g2  <- split(bio$exDB_ID, bio$Ortho_ID)
u   <- unique(bio[, .(Ortho_ID, or, p)]); l2 <- as.numeric(u$or>1 & u$p<0.05); names(l2) <- u$Ortho_ID
B   <- new("topGOdata", ontology="BP", allGenes=l2, geneSel=Multiple_threshold_test,
           annot=annFUN.gene2GO, gene2GO=g2, nodeSize=10)
rc  <- runTest(B, algorithm="classic",  statistic="fisher")
rw  <- runTest(B, algorithm="weight01", statistic="fisher")
tb  <- data.table(GenTable(B, classicFisher=rc, weight01=rw, orderBy="classicFisher", topNodes=2124))
fwrite(tb, "files/orthoDB/GO_enrich_protostomes_corrected.csv")

go <- head(tb, 12)
go[, p := as.numeric(gsub("<", "", classicFisher))]
go[, Significant := as.integer(Significant)][, Annotated := as.integer(Annotated)]
## DNA-repair-related terms among the corrected top 12
repair_terms <- tb[grepl("DNA repair|DNA damage|recombinational repair|double-strand break|nucleotide-excision|base-excision|mismatch repair|DNA metabolic", Term, ignore.case=TRUE)]$GO.ID
go[, repair := GO.ID %in% repair_terms]
suppressMessages(library(GO.db))
full <- AnnotationDbi::Term(GO.db::GOTERM[go$GO.ID])
go[, Term := factor(make.unique(unname(full)), levels = rev(make.unique(unname(full))))]
cat("\n--- corrected top 12 (classic Fisher) ---\n"); print(go[, .(GO.ID, Term, Significant, Annotated, p, repair)])

a <- ggplot(go, aes(-log10(p), Term, fill = repair)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%d/%d", Significant, Annotated)), hjust = -0.1, size = 1.6) +
  scale_fill_manual(values = c(`TRUE`="#d95f02", `FALSE`="grey70"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = expression(-log[10](italic(p))), y = NULL, title = "a")

## ---------- panel b: corrected Pagel / BH ----------
sig <- fread(file.path(out, "analysis/protostome_analysis/TableS5_corrected.csv"))
setorder(sig, P)
sig[, Name := factor(Name, levels = rev(Name))]
sig[, passBH := BH < 0.05]
cat("\npanel b genes:", nrow(sig), " passing BH:", sum(sig$passBH), "\n")

b <- ggplot(sig, aes(-log10(P), Name, fill = passBH)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = -log10(0.05), lty = 2, linewidth = 0.2) +
  scale_fill_manual(values = c(`TRUE`="#1b9e77", `FALSE`="grey75"),
                    labels = c(`TRUE`="BH-adjusted P < 0.05", `FALSE`="nominal P < 0.05 only"),
                    name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme(legend.position = "bottom", legend.key.size = unit(3,"mm"), legend.text = element_text(size=5)) +
  labs(x = expression(-log[10](italic(P)[Pagel])), y = NULL, title = "b")

fig <- (a / b) + plot_layout(heights = c(0.55, 1.15))
pdf(file.path(out, "revision/figures/FigS7_graphic.pdf"), width = 4.5, height = 5.4); print(fig); dev.off()
ggsave(file.path(out, "analysis/protostome_analysis/FigS_protostome_coevolution.png"), fig, width=4.5, height=5.4, dpi=600)
pdf(file.path(out, "analysis/protostome_analysis/FigS_protostome_coevolution.pdf"), width=4.5, height=5.4); print(fig); dev.off()
cat("\nwrote FigS7_graphic.pdf + analysis copies\n")
