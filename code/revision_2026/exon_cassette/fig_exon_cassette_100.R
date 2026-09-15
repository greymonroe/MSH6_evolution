## MSH6 domain architecture with exon boundaries — 76 species across the Tree of Life.

library(data.table)
library(ggplot2)
library(ape)
library(ggtree)
library(patchwork)

theme_set(theme_classic(base_size = 5))

hub <- "~/Dropbox/Research/MSH6_tudor/eLife_resubmission_2026/exon_cassette_analysis"

exons <- fread(file.path(hub, "exon_data_100.tsv"))
doms  <- fread(file.path(hub, "domain_data_100.tsv"))
summ  <- fread(file.path(hub, "species_summary_100.tsv"))

tree <- read.tree(file.path(hub, "species_tree_100.nwk"))
tree <- ladderize(tree, right = FALSE)

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

prot_len <- summ[, .(species, len = as.numeric(total_aa))]
prot_len[, species := factor(species, levels = sp_order)]

eb <- exons[exon_num > 1, .(species, pos = start_aa)]
eb[, species := factor(species, levels = sp_order)]

sp_labels <- gsub("_", " ", sp_order)
names(sp_labels) <- sp_order

p_tree <- ggtree(tree, linewidth = 0.2) +
  theme_tree2() +
  scale_x_continuous(labels = function(x) abs(x), name = "MYA") +
  theme(axis.text.x  = element_text(size = 3),
        axis.title.x = element_text(size = 4),
        plot.margin  = margin(2, 0, 2, 2, "pt"))

bar_h <- 0.4

p_doms <- ggplot() +
  geom_segment(data = prot_len,
               aes(x = 1, xend = len, y = species, yend = species),
               linewidth = 2, color = "grey85", lineend = "round") +
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
               linewidth = 0.15, color = "black") +
  scale_fill_manual(values = domain_cols, name = NULL,
                    guide = guide_legend(nrow = 1)) +
  scale_y_discrete(limits = sp_order, labels = sp_labels) +
  scale_x_continuous(breaks = seq(0, 1600, 200), name = "Amino acid position",
                     expand = expansion(mult = c(0.01, 0.02))) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 3, face = "italic"),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position  = "top",
    legend.key.size  = unit(0.2, "cm"),
    legend.text      = element_text(size = 4),
    legend.margin    = margin(0, 0, 0, 0),
    plot.margin      = margin(2, 4, 2, 0, "pt")
  )

fig <- p_tree + p_doms + plot_layout(widths = c(1, 4))

pdf(file.path(hub, "fig_exon_cassette_100.pdf"), width = 10, height = 11)
print(fig)
dev.off()

ggsave(file.path(hub, "fig_exon_cassette_100.png"), fig,
       width = 10, height = 11, dpi = 600)

cat("wrote fig_exon_cassette_100.{pdf,png}\n")
