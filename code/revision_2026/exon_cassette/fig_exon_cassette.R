## MSH6 domain architecture with exon boundaries across the Tree of Life.
## Data from: exon_data.tsv (NCBI API), domain_data.tsv (InterProScan), species_tree_41.nwk (TimeTree).
## Style matches Fig 1b: Tudor=green, PWWP=red, MutS=blue.

library(data.table)
library(ggplot2)
library(ape)
library(ggtree)
library(patchwork)

theme_set(theme_classic(base_size = 6))

hub <- "~/Dropbox/Research/MSH6_tudor/eLife_resubmission_2026/exon_cassette_analysis"

# ── Read data ────────────────────────────────────────────────────────────
exons  <- fread(file.path(hub, "exon_data.tsv"))
doms   <- fread(file.path(hub, "domain_data.tsv"))
summ   <- fread(file.path(hub, "species_summary.tsv"))

# ── Tree ─────────────────────────────────────────────────────────────────
tree <- read.tree(file.path(hub, "species_tree_41.nwk"))
tree <- ladderize(tree, right = FALSE)

gtree_plot <- ggtree(tree)
tip_data <- gtree_plot$data[gtree_plot$data$isTip, ]
tip_data <- tip_data[order(tip_data$y), ]
sp_order <- tip_data$label

# ── Clean domain names for display ───────────────────────────────────────
doms[, domain_label := gsub("_", " ", domain)]
doms[, domain_label := gsub("MutS (.*)", "MutS \\1", domain_label)]

# Collapse MutS subtypes for coloring: I,II,IV,V = same blue, III = purple
doms[, color_group := fcase(
  domain == "Tudor",   "Tudor",
  domain == "PWWP",    "PWWP",
  domain == "MutS_III","MutS III",
  grepl("MutS", domain), "MutS (other)"
)]

# ── Colors (matching Fig 1b) ─────────────────────────────────────────────
domain_cols <- c(
  "Tudor"        = "#4CAF50",
  "PWWP"         = "#E74C3C",
  "MutS (other)" = "#5B7DB1",
  "MutS III"     = "#7B68AE"
)

# ── Factor levels from tree tip order ────────────────────────────────────
exons[, species := factor(species, levels = sp_order)]
doms[,  species := factor(species, levels = sp_order)]
summ[,  species := factor(species, levels = sp_order)]

# ── Protein lengths ──────────────────────────────────────────────────────
prot_len <- summ[, .(species, len = as.numeric(total_aa))]
prot_len[, species := factor(species, levels = sp_order)]

# ── Exon boundaries (first aa of each exon, skip exon 1) ────────────────
eb <- exons[exon_num > 1, .(species, pos = start_aa)]
eb[, species := factor(species, levels = sp_order)]

# ── Pretty species labels ───────────────────────────────────────────────
sp_labels <- gsub("_", " ", sp_order)
names(sp_labels) <- sp_order

# ── Panel A: phylogenetic tree ───────────────────────────────────────────
p_tree <- ggtree(tree, size = 0.3) +
  theme_tree2() +
  scale_x_continuous(labels = function(x) abs(x), name = "MYA") +
  theme(axis.text.x = element_text(size = 4),
        axis.title.x = element_text(size = 5),
        plot.margin = margin(5, 0, 5, 5, "pt"))

# ── Panel B: domain architecture + exon boundaries ──────────────────────
bar_h <- 0.35

p_doms <- ggplot() +
  geom_segment(data = prot_len,
               aes(x = 1, xend = len, y = species, yend = species),
               linewidth = 1.5, color = "grey85", lineend = "round") +
  geom_rect(data = doms,
            aes(xmin = start_aa, xmax = end_aa,
                ymin = as.numeric(species) - bar_h,
                ymax = as.numeric(species) + bar_h,
                fill = color_group),
            color = NA) +
  geom_segment(data = eb,
               aes(x = pos, xend = pos,
                   y = as.numeric(species) - bar_h - 0.05,
                   yend = as.numeric(species) + bar_h + 0.05),
               linewidth = 0.15, color = "black") +
  scale_fill_manual(values = domain_cols, name = NULL,
                    guide = guide_legend(nrow = 1)) +
  scale_y_discrete(limits = sp_order, labels = sp_labels) +
  scale_x_continuous(breaks = seq(0, 1600, 200), name = "Amino acid position",
                     expand = expansion(mult = c(0.01, 0.02))) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 4, face = "italic"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.key.size = unit(0.25, "cm"),
    legend.text = element_text(size = 5),
    plot.margin = margin(5, 8, 5, 0, "pt")
  )

# ── Combine with patchwork ──────────────────────────────────────────────
fig <- p_tree + p_doms + plot_layout(widths = c(1, 2.5))

# ── Save ────────────────────────────────────────────────────────────────
pdf(file.path(hub, "fig_exon_cassette.pdf"), width = 7.5, height = 9.5)
print(fig)
dev.off()

ggsave(file.path(hub, "fig_exon_cassette.png"), fig,
       width = 7.5, height = 9.5, dpi = 600)

cat("wrote fig_exon_cassette.{pdf,png}\n")
