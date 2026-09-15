## MSH6 domain architecture, ALL placeable species, cassettes ALIGNED to the
## start of the PWWP / Tudor reader domain (x = 0 = reader start).
## Reader-less species anchored at N-terminus. Condensed to ~minimal row pitch
## with size-5 species labels, plus a colored major-clade strip.

library(data.table)
library(ggplot2)
library(ape)
library(ggtree)
library(patchwork)

theme_set(theme_classic(base_size = 5))

hub <- "~/Dropbox/Research/MSH6_tudor/eLife_resubmission_2026/analysis/exon_cassette_analysis"

PT_PER_TIP <- 4.8            # vertical points per species row
LAB_SIZE   <- 5             # species label font size (pt) — publication minimum

exons <- fread(file.path(hub, "exon_data_all.tsv"))
doms  <- fread(file.path(hub, "domain_data_all.tsv"))
summ  <- fread(file.path(hub, "species_summary_all.tsv"))

tree <- read.tree("~/repos/tol_reader_repair/files/time_tree_may4.nwk")

keep <- intersect(tree$tip.label, unique(summ$species))
coord_ck <- exons[species %in% keep, .(maxend = max(end_aa), minstart = min(start_aa)), by = species]
valid    <- coord_ck[maxend > 0 & maxend <= 3000 & minstart >= 0, species]
keep     <- intersect(keep, valid)

tree <- drop.tip(tree, setdiff(tree$tip.label, keep))
tree <- ladderize(tree, right = FALSE)

exons <- exons[species %in% keep]; doms <- doms[species %in% keep]; summ <- summ[species %in% keep]

# ---- alignment offset: start of PWWP/Tudor reader domain, else 0 ----
anchor <- doms[domain %in% c("PWWP", "Tudor"), .(offset = min(start_aa)), by = species]
off    <- setNames(rep(0, length(keep)), keep); off[anchor$species] <- anchor$offset
shift  <- function(sp, x) x - off[as.character(sp)]
cat("species:", length(keep), "| reader-anchored:", nrow(anchor), "\n")

# ---- major-clade assignment (specific -> broad; first match wins) ----
mrca_clade <- function(reps) {
  r <- reps[reps %in% tree$tip.label]
  if (length(r) < 2) return(character(0))
  extract.clade(tree, getMRCA(tree, r))$tip.label
}
clade_defs <- list(
  Vertebrata    = c("Homo_sapiens","Danio_rerio"),
  Insecta       = c("Drosophila_melanogaster","Apis_mellifera"),
  Arthropoda    = c("Drosophila_melanogaster","Limulus_polyphemus"),
  Nematoda      = c("Caenorhabditis_elegans","Caenorhabditis_briggsae"),
  Metazoa       = c("Homo_sapiens","Amphimedon_queenslandica"),
  Embryophyta   = c("Arabidopsis_thaliana","Selaginella_moellendorffii"),
  Viridiplantae = c("Arabidopsis_thaliana","Chlorella_vulgaris"),
  Fungi         = c("Schizosaccharomyces_pombe","Cryptococcus_neoformans")
)
clade_of <- setNames(rep("Other eukaryotes", length(keep)), keep)
assigned <- rep(FALSE, length(keep)); names(assigned) <- keep
for (cl in names(clade_defs)) {
  mem <- intersect(mrca_clade(clade_defs[[cl]]), keep)
  take <- mem[!assigned[mem]]
  clade_of[take] <- cl; assigned[take] <- TRUE
}
clade_levels <- c("Vertebrata","Metazoa","Insecta","Arthropoda","Nematoda",
                  "Fungi","Embryophyta","Viridiplantae","Other eukaryotes")
clade_cols <- c(Vertebrata="#1f78b4", Metazoa="#a6cee3", Insecta="#33a02c",
                Arthropoda="#b2df8a", Nematoda="#6a3d9a", Fungi="#ff7f00",
                Embryophyta="#e31a1c", Viridiplantae="#fb9a99",
                `Other eukaryotes`="#dddddd")

# ---- tip order ----
gtree_plot <- ggtree(tree)
tip_data <- gtree_plot$data[gtree_plot$data$isTip, ]
tip_data <- tip_data[order(tip_data$y), ]
sp_order <- tip_data$label

# ---- domain / exon / backbone geometry (shifted) ----
doms[, color_group := fcase(
  domain == "Tudor", "Tudor", domain == "PWWP", "PWWP",
  domain == "MutS_III", "MutS III", grepl("MutS", domain), "MutS (other)")]
domain_cols <- c("Tudor"="#4CAF50","PWWP"="#E74C3C","MutS (other)"="#5B7DB1","MutS III"="#7B68AE")

doms[, xmin := shift(species, start_aa)]; doms[, xmax := shift(species, end_aa)]
exons[, pos := shift(species, start_aa)]
prot_len <- exons[, .(len = max(end_aa)), by = species]
prot_len[, xstart := shift(species, 1)]; prot_len[, xend := shift(species, len)]
for (D in list(doms, exons, prot_len)) D[, species := factor(species, levels = sp_order)]
eb <- exons[exon_num > 1, .(species, pos)]
sp_labels <- gsub("_", " ", sp_order); names(sp_labels) <- sp_order

# clade strip data
strip <- data.table(species = factor(sp_order, levels = sp_order),
                    clade = factor(clade_of[sp_order], levels = clade_levels))
# clade name positions (vertical center of each clade band); label only sizeable clades
strip[, y := as.numeric(species)]
clade_lab <- strip[, .(y = mean(y), n = .N), by = clade][n >= 15 & clade != "Other eukaryotes"]
cat("\nclade band centers (tip index of", n_tip <- nrow(strip), "):\n"); print(clade_lab[order(y)])

n_tip <- length(sp_order)
fig_h <- n_tip * PT_PER_TIP / 72

# ---- panels ----
p_tree <- ggtree(tree, linewidth = 0.15) +
  theme_tree2() +
  scale_x_continuous(labels = function(x) abs(x), name = "MYA") +
  theme(axis.text.x = element_text(size = 4), axis.title.x = element_text(size = 5),
        plot.margin = margin(2, 0, 2, 2, "pt"))

p_clade <- ggplot(strip) +
  geom_rect(aes(xmin = 0, xmax = 1, ymin = y - 0.5, ymax = y + 0.5, fill = clade), color = NA) +
  geom_text(data = clade_lab, aes(x = 0.5, y = y, label = clade),
            angle = 90, size = 3.0, fontface = "bold", color = "black") +
  scale_fill_manual(values = clade_cols, name = "Clade", drop = FALSE,
                    guide = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  scale_y_continuous(limits = c(0.5, n_tip + 0.5), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() + theme(plot.margin = margin(2, 0, 2, 0, "pt"))

bar_h <- 0.42
p_doms <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.2, color = "grey40") +
  geom_segment(data = prot_len, aes(x = xstart, xend = xend, y = species, yend = species),
               linewidth = 0.9, color = "grey85", lineend = "round") +
  geom_rect(data = doms, aes(xmin = xmin, xmax = xmax,
                ymin = as.numeric(species) - bar_h, ymax = as.numeric(species) + bar_h,
                fill = color_group), color = NA) +
  geom_segment(data = eb, aes(x = pos, xend = pos,
                y = as.numeric(species) - bar_h - 0.03, yend = as.numeric(species) + bar_h + 0.03),
               linewidth = 0.1, color = "black") +
  scale_fill_manual(values = domain_cols, name = NULL, na.translate = FALSE,
                    guide = guide_legend(nrow = 1)) +
  scale_y_discrete(limits = sp_order, labels = sp_labels) +
  scale_x_continuous(breaks = seq(-400, 2200, 200),
                     name = "Position relative to PWWP/Tudor start (aa)  —  reader-less spp. anchored at N-terminus",
                     expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(xlim = c(-400, 2200)) +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_text(size = LAB_SIZE, face = "italic"),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top", legend.key.size = unit(0.25, "cm"),
        legend.text = element_text(size = 7), legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(2, 4, 2, 0, "pt"))

fig <- p_tree + p_clade + p_doms +
  plot_layout(widths = c(1, 0.22, 4), guides = "collect") &
  theme(legend.position = "top", legend.box = "horizontal")

pdf(file.path(hub, "fig_exon_cassette_aligned.pdf"), width = 12, height = fig_h)
print(fig); invisible(dev.off())
cat(sprintf("wrote fig_exon_cassette_aligned.pdf (%d species, %.1f in tall, %.1f pt/tip)\n",
            n_tip, fig_h, PT_PER_TIP))
print(table(clade_of))
