## MSH6 domain architecture with exon boundaries — ALL species with an exact TimeTree placement.
## Same layout as fig_exon_cassette_100.R but using the comprehensive *_all.tsv data,
## restricted to species that have an EXACT tip in the TimeTree (no genus grafting).

library(data.table)
library(ggplot2)
library(ape)
library(ggtree)
library(patchwork)

theme_set(theme_classic(base_size = 5))

hub <- "~/Dropbox/Research/MSH6_tudor/eLife_resubmission_2026/analysis/exon_cassette_analysis"

exons <- fread(file.path(hub, "exon_data_all.tsv"))
doms  <- fread(file.path(hub, "domain_data_all.tsv"))
summ  <- fread(file.path(hub, "species_summary_all.tsv"))

tree <- read.tree("~/repos/tol_reader_repair/files/time_tree_may4.nwk")

# Keep only species present (exactly) in both the data and the tree
keep <- intersect(tree$tip.label, unique(summ$species))

# The comprehensive *_all dataset contains ~120 species with corrupt (negative /
# absurd) exon coordinates; drop them so we only plot species with valid data.
coord_ck <- exons[species %in% keep, .(maxend = max(end_aa), minstart = min(start_aa)), by = species]
valid    <- coord_ck[maxend > 0 & maxend <= 3000 & minstart >= 0, species]
cat("Placeable species:", length(keep),
    "| clean:", length(valid),
    "| dropped (corrupt coords):", length(keep) - length(valid), "\n")
keep <- intersect(keep, valid)

tree <- drop.tip(tree, setdiff(tree$tip.label, keep))
tree <- ladderize(tree, right = FALSE)

exons <- exons[species %in% keep]
doms  <- doms[species %in% keep]
summ  <- summ[species %in% keep]

gtree_plot <- ggtree(tree)
tip_data <- gtree_plot$data[gtree_plot$data$isTip, ]
tip_data <- tip_data[order(tip_data$y), ]
sp_order <- tip_data$label

doms[, color_group := fcase(
  domain == "Tudor",    "Tudor",
  domain == "PWWP",     "PWWP",
  domain == "MutS_III", "MutS III",
  grepl("MutS", domain), "MutS (other)"
)]

domain_cols <- c(
  "Tudor"        = "#4CAF50",
  "PWWP"         = "#E74C3C",
  "MutS (other)" = "#5B7DB1",
  "MutS III"     = "#7B68AE"
)

exons[, species := factor(species, levels = sp_order)]
doms[,  species := factor(species, levels = sp_order)]
summ[,  species := factor(species, levels = sp_order)]

# Backbone length from exon coordinates (consistent with plotted domains/exons,
# and immune to the corrupt total_aa column in species_summary_all.tsv)
prot_len <- exons[, .(len = max(end_aa)), by = species]
prot_len[, species := factor(species, levels = sp_order)]

eb <- exons[exon_num > 1, .(species, pos = start_aa)]
eb[, species := factor(species, levels = sp_order)]

sp_labels <- gsub("_", " ", sp_order)
names(sp_labels) <- sp_order

n_tip <- length(sp_order)
fig_h <- max(11, n_tip * 0.085)   # ~0.085 in per tip

p_tree <- ggtree(tree, linewidth = 0.15) +
  theme_tree2() +
  scale_x_continuous(labels = function(x) abs(x), name = "MYA") +
  theme(axis.text.x  = element_text(size = 3),
        axis.title.x = element_text(size = 4),
        plot.margin  = margin(2, 0, 2, 2, "pt"))

bar_h <- 0.4

p_doms <- ggplot() +
  geom_segment(data = prot_len,
               aes(x = 1, xend = len, y = species, yend = species),
               linewidth = 1.2, color = "grey85", lineend = "round") +
  geom_rect(data = doms,
            aes(xmin = start_aa, xmax = end_aa,
                ymin = as.numeric(species) - bar_h,
                ymax = as.numeric(species) + bar_h,
                fill = color_group),
            color = NA) +
  geom_segment(data = eb,
               aes(x = pos, xend = pos,
                   y = as.numeric(species) - bar_h - 0.03,
                   yend = as.numeric(species) + bar_h + 0.03),
               linewidth = 0.1, color = "black") +
  scale_fill_manual(values = domain_cols, name = NULL, na.translate = FALSE,
                    guide = guide_legend(nrow = 1)) +
  scale_y_discrete(limits = sp_order, labels = sp_labels) +
  scale_x_continuous(breaks = seq(0, 2400, 200), name = "Amino acid position",
                     expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(xlim = c(0, 2500)) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 2, face = "italic"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position  = "top",
    legend.key.size  = unit(0.2, "cm"),
    legend.text      = element_text(size = 5),
    legend.margin    = margin(0, 0, 0, 0),
    plot.margin      = margin(2, 4, 2, 0, "pt")
  )

fig <- p_tree + p_doms + plot_layout(widths = c(1, 4))

pdf(file.path(hub, "fig_exon_cassette_all.pdf"), width = 11, height = fig_h)
print(fig)
invisible(dev.off())

cat(sprintf("wrote fig_exon_cassette_all.pdf (%d species, %.0f in tall)\n", n_tip, fig_h))
