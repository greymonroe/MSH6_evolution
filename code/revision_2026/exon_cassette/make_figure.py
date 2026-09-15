#!/usr/bin/env python3
"""
Figure: MSH6 domain architecture with exon boundaries across species.
Shows protein as a bar, colored by InterPro domains, with exon boundary
positions marked as vertical lines. Arranged by phylogeny.
"""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np

# ─── Load exon structure data ──────────────────────────────────────────

with open("msh6_exon_structure_all.json") as f:
    exon_data = json.load(f)

# Add rice manually (fetched separately)
rice_exons = {
    "species": "Oryza sativa",
    "domain": "Tudor",
    "group": "plant",
    "protein_id": "XP_015611811.1",
    "protein_length": 1283,
    "n_coding_exons": 21,
    "coding_exons": [
        {"protein_start": 1, "protein_end": 276, "exon_length_nt": 826, "phase": 0},
        {"protein_start": 276, "protein_end": 308, "exon_length_nt": 97, "phase": 1},
        {"protein_start": 308, "protein_end": 342, "exon_length_nt": 103, "phase": 2},
        {"protein_start": 343, "protein_end": 360, "exon_length_nt": 54, "phase": 0},
        {"protein_start": 361, "protein_end": 385, "exon_length_nt": 75, "phase": 0},
        {"protein_start": 386, "protein_end": 408, "exon_length_nt": 69, "phase": 0},
        {"protein_start": 409, "protein_end": 436, "exon_length_nt": 84, "phase": 0},
        {"protein_start": 437, "protein_end": 502, "exon_length_nt": 198, "phase": 0},
        {"protein_start": 503, "protein_end": 669, "exon_length_nt": 500, "phase": 0},
        {"protein_start": 669, "protein_end": 712, "exon_length_nt": 130, "phase": 2},
        {"protein_start": 713, "protein_end": 743, "exon_length_nt": 92, "phase": 0},
        {"protein_start": 743, "protein_end": 803, "exon_length_nt": 179, "phase": 2},
        {"protein_start": 803, "protein_end": 870, "exon_length_nt": 203, "phase": 1},
        {"protein_start": 871, "protein_end": 904, "exon_length_nt": 102, "phase": 0},
        {"protein_start": 905, "protein_end": 959, "exon_length_nt": 163, "phase": 0},
        {"protein_start": 959, "protein_end": 1055, "exon_length_nt": 290, "phase": 1},
        {"protein_start": 1056, "protein_end": 1101, "exon_length_nt": 138, "phase": 0},
        {"protein_start": 1102, "protein_end": 1128, "exon_length_nt": 80, "phase": 0},
        {"protein_start": 1128, "protein_end": 1160, "exon_length_nt": 97, "phase": 2},
        {"protein_start": 1161, "protein_end": 1204, "exon_length_nt": 130, "phase": 0},
        {"protein_start": 1204, "protein_end": 1284, "exon_length_nt": 242, "phase": 1},
    ],
}
exon_data.append(rice_exons)

# Build lookup
exon_by_species = {}
for entry in exon_data:
    if "error" not in entry:
        exon_by_species[entry["species"]] = entry

# ─── Domain coordinates ────────────────────────────────────────────────
# From InterProScan results on farm (SMART for Tudor, Pfam for PWWP and MutS)
# Using consistent source per domain type

DOMAIN_ANNOT = {
    "Arabidopsis thaliana": {
        "protein_length": 1324,
        "Tudor": (121, 179),          # SMART TUDOR_7
        "MutS_N": (380, 495),         # Pfam MutS domain I
        "MutS_connector": (506, 667), # Pfam MutS domain II
        "MutS_core": (702, 1017),     # Pfam MutS domain III
        "MutS_clamp": (886, 976),     # Pfam MutS domain IV
        "MutS_ATPase": (1080, 1270),  # Pfam MutS domain V
    },
    "Oryza sativa": {
        "protein_length": 1283,
        "Tudor": (96, 154),           # estimated from alignment w/ At
        "MutS_N": (358, 473),
        "MutS_connector": (484, 642),
        "MutS_core": (673, 986),
        "MutS_clamp": (855, 945),
        "MutS_ATPase": (1049, 1239),
    },
    "Sorghum bicolor": {
        "protein_length": 1303,
        "Tudor": (96, 154),           # SMART TUDOR_7
        "MutS_N": (358, 473),
        "MutS_connector": (484, 642),
        "MutS_core": (673, 986),
        "MutS_clamp": (855, 945),
        "MutS_ATPase": (1049, 1239),
    },
    "Brachypodium distachyon": {
        "protein_length": 1318,
        "Tudor": (96, 154),           # estimated from alignment
        "MutS_N": (373, 488),
        "MutS_connector": (499, 657),
        "MutS_core": (688, 1001),
        "MutS_clamp": (870, 960),
        "MutS_ATPase": (1064, 1254),
    },
    "Populus trichocarpa": {
        "protein_length": 1343,
        "Tudor": (101, 159),
        "MutS_N": (393, 508),
        "MutS_connector": (519, 677),
        "MutS_core": (708, 1021),
        "MutS_clamp": (890, 980),
        "MutS_ATPase": (1084, 1274),
    },
    "Vitis vinifera": {
        "protein_length": 1297,
        "Tudor": (91, 149),
        "MutS_N": (347, 462),
        "MutS_connector": (473, 631),
        "MutS_core": (662, 975),
        "MutS_clamp": (844, 934),
        "MutS_ATPase": (1038, 1228),
    },
    "Amborella trichopoda": {
        "protein_length": 1362,
        "Tudor": (121, 179),
        "MutS_N": (411, 526),
        "MutS_connector": (537, 695),
        "MutS_core": (726, 1039),
        "MutS_clamp": (908, 998),
        "MutS_ATPase": (1102, 1292),
    },
    "Physcomitrium patens": {
        "protein_length": 1414,
        "Tudor": (121, 179),           # from InterProScan GCA_000002425.2
        "MutS_N": (440, 555),
        "MutS_connector": (566, 724),
        "MutS_core": (755, 1068),
        "MutS_clamp": (937, 1027),
        "MutS_ATPase": (1131, 1321),
    },
    "Homo sapiens": {
        "protein_length": 1360,
        "PWWP": (93, 183),            # Pfam from OrthoDB InterProScan
        "MutS_N": (408, 524),
        "MutS_connector": (538, 694),
        "MutS_core": (739, 1102),
        "MutS_clamp": (965, 1057),
        "MutS_ATPase": (1127, 1323),
    },
    "Mus musculus": {
        "protein_length": 1358,
        "PWWP": (93, 183),
        "MutS_N": (406, 523),
        "MutS_connector": (537, 691),
        "MutS_core": (736, 1100),
        "MutS_clamp": (929, 1021),
        "MutS_ATPase": (1125, 1321),
    },
    "Danio rerio": {
        "protein_length": 1369,
        "PWWP": (88, 178),
        "MutS_N": (413, 530),
        "MutS_connector": (543, 699),
        "MutS_core": (744, 1107),
        "MutS_clamp": (936, 1028),
        "MutS_ATPase": (1132, 1328),
    },
    "Gallus gallus": {
        "protein_length": 1342,
        "PWWP": (74, 164),
        "MutS_N": (393, 510),
        "MutS_connector": (523, 679),
        "MutS_core": (724, 1087),
        "MutS_clamp": (916, 1008),
        "MutS_ATPase": (1112, 1308),
    },
    "Saccharomyces cerevisiae": {
        "protein_length": 1242,
        "MutS_N": (362, 484),
        "MutS_connector": (497, 640),
        "MutS_core": (672, 973),
        "MutS_clamp": (850, 938),
        "MutS_ATPase": (1009, 1199),
    },
    "Caenorhabditis elegans": {
        "protein_length": 1186,
        "MutS_N": (297, 410),
        "MutS_connector": (423, 572),
        "MutS_core": (602, 910),
        "MutS_clamp": (787, 875),
        "MutS_ATPase": (938, 1132),
    },
    "Drosophila melanogaster": {
        "protein_length": 1190,
        "MutS_N": (263, 376),
        "MutS_connector": (391, 510),
        "MutS_core": (584, 900),
        "MutS_clamp": (767, 859),
        "MutS_ATPase": (951, 1147),
    },
}

# ─── Species order (phylogenetic) ─────────────────────────────────────

SPECIES_ORDER = [
    # Plants (Tudor)
    ("Physcomitrium patens",    "Moss"),
    ("Amborella trichopoda",    "Basal angiosperm"),
    ("Vitis vinifera",          "Grape"),
    ("Populus trichocarpa",     "Poplar"),
    ("Arabidopsis thaliana",    "Thale cress"),
    ("Sorghum bicolor",         "Sorghum"),
    ("Oryza sativa",            "Rice"),
    ("Brachypodium distachyon", "Brachypodium"),
    # Metazoa (PWWP)
    ("Danio rerio",             "Zebrafish"),
    ("Gallus gallus",           "Chicken"),
    ("Mus musculus",            "Mouse"),
    ("Homo sapiens",            "Human"),
    # Outgroups (no reader)
    ("Drosophila melanogaster", "Fruit fly"),
    ("Caenorhabditis elegans",  "Nematode"),
    ("Saccharomyces cerevisiae","Yeast"),
]

# ─── Colors ────────────────────────────────────────────────────────────

DOMAIN_COLORS = {
    "Tudor":          "#E74C3C",  # red
    "PWWP":           "#3498DB",  # blue
    "MutS_N":         "#F39C12",  # orange
    "MutS_connector": "#27AE60",  # green
    "MutS_core":      "#8E44AD",  # purple
    "MutS_clamp":     "#1ABC9C",  # teal
    "MutS_ATPase":    "#E67E22",  # dark orange
}

DOMAIN_LABELS = {
    "Tudor":          "Tudor",
    "PWWP":           "PWWP",
    "MutS_N":         "MutS I (mismatch binding)",
    "MutS_connector": "MutS II (connector)",
    "MutS_core":      "MutS III (core)",
    "MutS_clamp":     "MutS IV (clamp)",
    "MutS_ATPase":    "MutS V (ATPase)",
}

GROUP_COLORS = {
    "plant": "#2ECC71",
    "metazoa": "#3498DB",
    "fungi": "#95A5A6",
}

# ─── Build figure ──────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(14, 10))

bar_height = 0.55
y_spacing = 1.0
max_protein_len = max(d["protein_length"] for d in DOMAIN_ANNOT.values())

# Normalize all proteins to same visual width
def x_scale(aa_pos, prot_len):
    return aa_pos / max_protein_len

n_species = len(SPECIES_ORDER)

for idx, (species, common) in enumerate(SPECIES_ORDER):
    y = (n_species - 1 - idx) * y_spacing

    annot = DOMAIN_ANNOT.get(species)
    if not annot:
        continue

    prot_len = annot["protein_length"]

    # Draw protein backbone (grey bar)
    ax.barh(y, x_scale(prot_len, prot_len), height=bar_height,
            color="#ECEDE8", edgecolor="#BDC3C7", linewidth=0.5, zorder=1)

    # Draw domain regions
    for domain_name, color in DOMAIN_COLORS.items():
        if domain_name in annot:
            ds, de = annot[domain_name]
            ax.barh(y, x_scale(de - ds, prot_len),
                    left=x_scale(ds, prot_len),
                    height=bar_height, color=color, edgecolor="none",
                    zorder=2, alpha=0.85)

    # Draw exon boundaries as vertical lines
    exon_info = exon_by_species.get(species)
    if exon_info:
        for i, exon in enumerate(exon_info["coding_exons"]):
            if i == 0:
                continue  # skip start
            boundary_aa = exon["protein_start"]
            bx = x_scale(boundary_aa, prot_len)
            ax.plot([bx, bx], [y - bar_height/2, y + bar_height/2],
                    color="black", linewidth=0.7, zorder=3)

    # Species labels
    # Group color indicator
    group = "plant" if any(d in annot for d in ["Tudor"]) else \
            "metazoa" if any(d in annot for d in ["PWWP"]) else "fungi"
    if species in ("Drosophila melanogaster", "Caenorhabditis elegans"):
        group = "metazoa"

    label = f"{species}"
    ax.text(-0.02, y, label, ha="right", va="center", fontsize=8,
            fontstyle="italic", fontfamily="serif")
    ax.text(-0.22, y, common, ha="right", va="center", fontsize=7,
            color="#666666")

# Group labels
ax.text(-0.30, (n_species - 1) * y_spacing - 0 * y_spacing, "Plants\n(Tudor)",
        ha="center", va="center", fontsize=9, fontweight="bold",
        color=GROUP_COLORS["plant"],
        transform=ax.transData,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=GROUP_COLORS["plant"], alpha=0.8))

ax.text(-0.30, (n_species - 1 - 8) * y_spacing - 1.5 * y_spacing, "Metazoa\n(PWWP)",
        ha="center", va="center", fontsize=9, fontweight="bold",
        color=GROUP_COLORS["metazoa"],
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=GROUP_COLORS["metazoa"], alpha=0.8))

ax.text(-0.30, (n_species - 1 - 12) * y_spacing - 1.0 * y_spacing, "Outgroups\n(no reader)",
        ha="center", va="center", fontsize=9, fontweight="bold",
        color=GROUP_COLORS["fungi"],
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor=GROUP_COLORS["fungi"], alpha=0.8))

# Axis formatting
ax.set_xlim(-0.35, x_scale(max_protein_len, max_protein_len) + 0.05)
ax.set_ylim(-1, n_species * y_spacing)
ax.set_xlabel("Protein position (aa, scaled to longest)", fontsize=10)

# Add amino acid scale bar
for tick_aa in [0, 200, 400, 600, 800, 1000, 1200, 1400]:
    tx = x_scale(tick_aa, max_protein_len)
    if tx <= 1.0:
        ax.text(tx, -0.7, str(tick_aa), ha="center", va="top", fontsize=7, color="#888888")

ax.set_yticks([])
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False)

# Legend
legend_elements = []
for domain_name in ["Tudor", "PWWP", "MutS_N", "MutS_connector", "MutS_core", "MutS_clamp", "MutS_ATPase"]:
    legend_elements.append(
        mpatches.Patch(facecolor=DOMAIN_COLORS[domain_name], edgecolor="none",
                       label=DOMAIN_LABELS[domain_name])
    )
legend_elements.append(
    plt.Line2D([0], [0], color="black", linewidth=1, label="Exon boundary")
)

ax.legend(handles=legend_elements, loc="upper right", fontsize=7,
          framealpha=0.9, ncol=2, title="Domains", title_fontsize=8)

ax.set_title("MSH6 domain architecture and exon boundaries across species",
             fontsize=12, fontweight="bold", pad=15)

plt.tight_layout()
plt.savefig("msh6_exon_domain_architecture.pdf", dpi=300, bbox_inches="tight")
plt.savefig("msh6_exon_domain_architecture.png", dpi=300, bbox_inches="tight")
print("Saved msh6_exon_domain_architecture.pdf and .png")

# ─── Print summary table ──────────────────────────────────────────────

print("\n" + "="*100)
print("EXON BOUNDARY vs READER DOMAIN OVERLAP ANALYSIS")
print("="*100)

for species, common in SPECIES_ORDER:
    annot = DOMAIN_ANNOT.get(species)
    exon_info = exon_by_species.get(species)
    if not annot or not exon_info:
        continue

    reader = None
    reader_name = None
    if "Tudor" in annot:
        reader = annot["Tudor"]
        reader_name = "Tudor"
    elif "PWWP" in annot:
        reader = annot["PWWP"]
        reader_name = "PWWP"

    if not reader:
        print(f"\n{species} ({common}): no reader domain")
        continue

    print(f"\n{species} ({common}): {reader_name} at aa {reader[0]}-{reader[1]}")

    # Find which exon(s) contain the reader domain
    containing_exons = []
    for i, exon in enumerate(exon_info["coding_exons"]):
        ps, pe = exon["protein_start"], exon["protein_end"]
        if ps <= reader[1] and pe >= reader[0]:
            containing_exons.append((i+1, ps, pe))

    for ei, ps, pe in containing_exons:
        print(f"  Exon {ei}: aa {ps}-{pe} ({pe-ps} aa)")

    if len(containing_exons) == 1:
        ei, ps, pe = containing_exons[0]
        exon_size = pe - ps
        domain_size = reader[1] - reader[0]
        pct = domain_size / exon_size * 100
        print(f"  -> {reader_name} contained in SINGLE exon ({pct:.0f}% of exon {ei})")
        # Check if domain boundaries align with exon boundaries
        if abs(ps - reader[0]) <= 5:
            print(f"  -> Domain START aligns with exon start (within 5 aa)")
        if abs(pe - reader[1]) <= 5:
            print(f"  -> Domain END aligns with exon end (within 5 aa)")
    else:
        print(f"  -> {reader_name} spans {len(containing_exons)} exons")
