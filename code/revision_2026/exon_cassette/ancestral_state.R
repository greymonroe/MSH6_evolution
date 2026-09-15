library(ape)
library(phytools)
library(data.table)

hub <- "~/Dropbox/Research/MSH6_tudor/eLife_resubmission_2026/exon_cassette_analysis"

tree <- read.tree("~/repos/tol_reader_repair/files/time_tree_may4.nwk")

traits <- fread(file.path(hub, "reader_exon_traits.tsv"))

# Handle tree-tip name mismatches (strain suffixes etc.)
tree_tips <- tree$tip.label
trait_tips <- traits$tree_tip

# Prune tree to species with trait data
keep <- intersect(tree_tips, trait_tips)
cat("Tips in both tree and data:", length(keep), "\n")

drop <- setdiff(tree_tips, keep)
tree_pruned <- drop.tip(tree, drop)
cat("Pruned tree tips:", Ntip(tree_pruned), "\n")

traits_matched <- traits[tree_tip %in% tree_pruned$tip.label]
trait_vec <- setNames(traits_matched$trait, traits_matched$tree_tip)
trait_vec <- trait_vec[tree_pruned$tip.label]

# ============================================================
# ANALYSIS 1: PWWP exon architecture
# Focus on PWWP species: is 2-exon the ancestral state?
# Recode: PWWP_2_exon, PWWP_single_exon, PWWP_3plus_exon, other
# ============================================================

# For PWWP-focused analysis, recode non-PWWP as "no_PWWP"
pwwp_trait <- ifelse(grepl("^PWWP", trait_vec), trait_vec, "no_PWWP")
pwwp_trait <- gsub("PWWP_unknown", "no_PWWP", pwwp_trait)
pwwp_trait <- factor(pwwp_trait)
names(pwwp_trait) <- names(trait_vec)

cat("\nPWWP trait distribution:\n")
print(table(pwwp_trait))

# Fit equal-rates model
cat("\n=== Fitting ER model for PWWP ===\n")
pwwp_er <- fitMk(tree_pruned, pwwp_trait, model = "ER")
cat("ER log-likelihood:", pwwp_er$logLik, "\n")

# Fit all-rates-different model
cat("\n=== Fitting ARD model for PWWP ===\n")
pwwp_ard <- fitMk(tree_pruned, pwwp_trait, model = "ARD")
cat("ARD log-likelihood:", pwwp_ard$logLik, "\n")

# AIC comparison
aic_er <- -2 * pwwp_er$logLik + 2 * length(pwwp_er$rates)
aic_ard <- -2 * pwwp_ard$logLik + 2 * length(pwwp_ard$rates)
cat("\nAIC ER:", aic_er, "  AIC ARD:", aic_ard, "\n")
cat("Best model:", ifelse(aic_ard < aic_er, "ARD", "ER"), "\n")

best_model <- if (aic_ard < aic_er) pwwp_ard else pwwp_er

# Ancestral state reconstruction using ape::ace (more stable API)
cat("\n=== Ancestral state probabilities at key nodes ===\n")
asr <- ace(pwwp_trait, tree_pruned, type = "discrete", model = "ARD")

# Root node
root_node <- Ntip(tree_pruned) + 1
cat("\nRoot node probabilities:\n")
root_probs <- asr$lik.anc[1, ]
print(round(root_probs, 4))

# Find MRCA of key clades for node probabilities
# We need to identify internal nodes at key positions

# Helper: find MRCA of two tips
find_mrca <- function(tip1, tip2) {
    t1_matches <- grep(tip1, tree_pruned$tip.label, value = TRUE)
    t2_matches <- grep(tip2, tree_pruned$tip.label, value = TRUE)
    if (length(t1_matches) == 0 || length(t2_matches) == 0) return(NA)
    getMRCA(tree_pruned, c(t1_matches[1], t2_matches[1]))
}

# Key ancestral nodes
nodes_to_check <- list(
    "MRCA Bilateria (Human+Drosophila)" = find_mrca("Homo_sapiens", "Drosophila"),
    "MRCA Deuterostomia (Human+Strongylocentrotus)" = find_mrca("Homo_sapiens", "Strongylocentrotus"),
    "MRCA Vertebrata (Human+Petromyzon)" = find_mrca("Homo_sapiens", "Petromyzon"),
    "MRCA Metazoa (Human+Amphimedon)" = find_mrca("Homo_sapiens", "Amphimedon"),
    "MRCA Ecdysozoa (Drosophila+C.elegans)" = find_mrca("Drosophila", "Caenorhabditis"),
    "MRCA Viridiplantae (Arabidopsis+Chlamydomonas)" = find_mrca("Arabidopsis", "Chlamydomonas")
)

for (label in names(nodes_to_check)) {
    node <- nodes_to_check[[label]]
    if (is.na(node)) {
        cat(label, ": node not found\n")
        next
    }
    node_idx <- node - Ntip(tree_pruned)
    if (node_idx < 1 || node_idx > nrow(asr$lik.anc)) {
        cat(label, ": index out of range\n")
        next
    }
    probs <- asr$lik.anc[node_idx, ]
    cat("\n", label, " (node ", node, "):\n", sep = "")
    for (state in names(probs)) {
        if (probs[state] > 0.01) {
            cat("  ", state, ": ", round(probs[state], 3), "\n", sep = "")
        }
    }
}

# ============================================================
# ANALYSIS 2: Simplified binary — does PWWP have internal exon break?
# States: "break" (2+ exons) vs "no_break" (single exon or no PWWP)
# This tests: was the ancestral PWWP a multi-exon or single-exon insertion?
# ============================================================
cat("\n\n=== Simplified: PWWP exon break (yes/no) among PWWP species only ===\n")

# Prune to only PWWP species
pwwp_only <- trait_vec[grepl("^PWWP_[23s]", trait_vec)]
pwwp_binary <- ifelse(grepl("2_exon|3plus", pwwp_only), "break", "no_break")
pwwp_binary <- factor(pwwp_binary)
names(pwwp_binary) <- names(pwwp_only)

tree_pwwp <- drop.tip(tree_pruned, setdiff(tree_pruned$tip.label, names(pwwp_binary)))
pwwp_binary <- pwwp_binary[tree_pwwp$tip.label]

cat("PWWP-only species:", Ntip(tree_pwwp), "\n")
cat("Trait counts:\n")
print(table(pwwp_binary))

if (Ntip(tree_pwwp) > 10 && length(levels(pwwp_binary)) == 2) {
    asr_binary <- ace(pwwp_binary, tree_pwwp, type = "discrete", model = "ARD")
    root_binary <- asr_binary$lik.anc[1, ]
    cat("\nRoot of PWWP clade:\n")
    print(round(root_binary, 4))
    cat("\nTransition rates (Q matrix):\n")
    print(round(asr_binary$rates, 6))
    cat("\nLog-likelihood:", asr_binary$loglik, "\n")
}

cat("\n=== Done ===\n")
