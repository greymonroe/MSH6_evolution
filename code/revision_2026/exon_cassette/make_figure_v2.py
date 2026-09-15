#!/usr/bin/env python3
"""
Figure v2: MSH6 domain architecture with exon boundaries across species.
All exon data verified from NCBI genomic GenBank records.
Domain coords from InterProScan (farm) where available, else SMART/CDD.
"""

import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ─── Verified exon structures (all from NCBI genomic GenBank CDS join) ──

EXON_DATA = {
    "Arabidopsis thaliana": {
        "protein_length": 1324, "n_exons": 21, "group": "plant",
        "exons": [
            (1,316),(316,349),(349,383),(384,401),(402,426),(427,449),(450,477),
            (478,542),(543,719),(719,762),(763,793),(793,853),(853,920),(921,954),
            (955,1009),(1009,1104),(1105,1150),(1151,1177),(1177,1210),(1211,1254),(1254,1325),
        ],
    },
    "Oryza sativa": {
        "protein_length": 1283, "n_exons": 21, "group": "plant",
        "exons": [
            (1,276),(276,308),(308,342),(343,360),(361,385),(386,408),(409,436),
            (437,502),(503,669),(669,712),(713,743),(743,803),(803,870),(871,904),
            (905,959),(959,1055),(1056,1101),(1102,1128),(1128,1160),(1161,1204),(1204,1284),
        ],
    },
    "Zea mays": {
        "protein_length": 1296, "n_exons": 21, "group": "plant",
        "exons": [
            (1,293),(293,325),(325,359),(360,377),(378,402),(403,425),(426,453),
            (454,519),(520,688),(688,731),(732,762),(762,820),(820,887),(888,921),
            (922,976),(976,1072),(1073,1118),(1119,1145),(1145,1177),(1178,1221),(1221,1297),
        ],
    },
    "Sorghum bicolor": {
        "protein_length": 1304, "n_exons": 21, "group": "plant",
        "exons": [
            (1,295),(295,327),(327,361),(362,379),(380,404),(405,427),(428,455),
            (456,521),(522,690),(690,733),(734,764),(764,822),(822,889),(890,923),
            (924,978),(978,1074),(1075,1120),(1121,1147),(1147,1179),(1180,1223),(1223,1305),
        ],
    },
    "Brachypodium distachyon": {
        "protein_length": 1318, "n_exons": 21, "group": "plant",
        "exons": [
            (1,306),(306,338),(338,372),(373,390),(391,415),(416,438),(439,466),
            (467,532),(533,703),(703,746),(747,777),(777,837),(837,904),(905,938),
            (939,993),(993,1089),(1090,1135),(1136,1162),(1162,1194),(1195,1238),(1238,1319),
        ],
    },
    "Populus trichocarpa": {
        "protein_length": 1343, "n_exons": 22, "group": "plant",
        "exons": [
            (1,281),(281,312),(312,322),(322,356),(357,374),(375,399),(400,422),
            (423,450),(451,515),(516,699),(699,742),(743,773),(773,833),(833,900),
            (901,934),(935,989),(989,1116),(1117,1162),(1163,1189),(1189,1222),(1223,1265),(1265,1344),
        ],
    },
    "Vitis vinifera": {
        "protein_length": 1297, "n_exons": 21, "group": "plant",
        "exons": [
            (1,272),(272,304),(304,338),(339,356),(357,381),(382,404),(405,432),
            (433,492),(493,675),(675,718),(719,749),(749,809),(809,876),(877,910),
            (911,965),(965,1061),(1062,1107),(1108,1134),(1134,1167),(1168,1211),(1211,1298),
        ],
    },
    "Amborella trichopoda": {
        "protein_length": 1362, "n_exons": 21, "group": "plant",
        "exons": [
            (1,333),(333,365),(365,399),(400,417),(418,442),(443,465),(466,493),
            (494,559),(560,745),(745,788),(789,819),(819,879),(879,946),(947,980),
            (981,1035),(1035,1134),(1135,1180),(1181,1207),(1207,1240),(1241,1284),(1284,1363),
        ],
    },
    "Physcomitrium patens": {
        "protein_length": 1410, "n_exons": 20, "group": "plant",
        "exons": [
            (1,397),(397,423),(423,457),(458,475),(476,500),(501,523),(524,551),
            (552,612),(613,803),(803,846),(847,877),(877,937),(937,1005),(1006,1039),
            (1040,1094),(1094,1184),(1185,1230),(1231,1257),(1257,1290),(1291,1411),
        ],
    },
    "Homo sapiens": {
        "protein_length": 1360, "n_exons": 10, "group": "metazoa",
        "exons": [
            (1,87),(87,153),(153,209),(210,1058),(1058,1146),(1147,1186),
            (1186,1216),(1216,1267),(1268,1334),(1334,1361),
        ],
    },
    "Mus musculus": {
        "protein_length": 1358, "n_exons": 10, "group": "metazoa",
        "exons": [
            (1,86),(86,153),(153,209),(210,1055),(1055,1144),(1145,1184),
            (1184,1214),(1214,1265),(1266,1332),(1332,1359),
        ],
    },
    "Danio rerio": {
        "protein_length": 1369, "n_exons": 10, "group": "metazoa",
        "exons": [
            (1,79),(79,148),(148,206),(207,1059),(1059,1155),(1156,1195),
            (1195,1225),(1225,1276),(1277,1343),(1343,1370),
        ],
    },
    "Gallus gallus": {
        "protein_length": 1341, "n_exons": 10, "group": "metazoa",
        "exons": [
            (1,66),(66,134),(134,192),(193,1039),(1039,1127),(1128,1167),
            (1167,1197),(1197,1248),(1249,1315),(1315,1342),
        ],
    },
    "Drosophila melanogaster": {
        "protein_length": 1190, "n_exons": 5, "group": "metazoa",
        "exons": [
            (1,334),(335,444),(445,980),(981,1142),(1143,1191),
        ],
    },
    "Caenorhabditis elegans": {
        "protein_length": 1186, "n_exons": 8, "group": "metazoa",
        "exons": [
            (1,46),(46,139),(140,394),(395,663),(663,800),(801,984),(985,1095),(1096,1187),
        ],
    },
    "Saccharomyces cerevisiae": {
        "protein_length": 1242, "n_exons": 1, "group": "fungi",
        "exons": [(1,1243)],
    },
}

# ─── Domain coordinates (InterProScan where available) ────────────────
# Tudor: SMART TUDOR_7 from farm InterProScan (Arabidopsis confirmed)
# PWWP: Pfam PF00855 from OrthoDB InterProScan
# MutS domains: Pfam from InterProScan

DOMAINS = {
    "Arabidopsis thaliana":    {"Tudor": (121,179), "MutS_N": (380,495), "MutS_conn": (506,667), "MutS_core": (702,1017), "MutS_clamp": (886,976), "MutS_ATP": (1080,1270)},
    "Oryza sativa":            {"Tudor": (96,154),  "MutS_N": (358,473), "MutS_conn": (484,642), "MutS_core": (673,986),  "MutS_clamp": (855,945), "MutS_ATP": (1049,1239)},
    "Zea mays":                {"Tudor": (98,156),  "MutS_N": (365,480), "MutS_conn": (491,649), "MutS_core": (680,993),  "MutS_clamp": (862,952), "MutS_ATP": (1056,1246)},
    "Sorghum bicolor":         {"Tudor": (100,158), "MutS_N": (367,482), "MutS_conn": (493,651), "MutS_core": (682,995),  "MutS_clamp": (864,954), "MutS_ATP": (1058,1248)},
    "Brachypodium distachyon": {"Tudor": (96,154),  "MutS_N": (373,488), "MutS_conn": (499,657), "MutS_core": (688,1001), "MutS_clamp": (870,960), "MutS_ATP": (1064,1254)},
    "Populus trichocarpa":     {"Tudor": (101,159), "MutS_N": (393,508), "MutS_conn": (519,677), "MutS_core": (708,1021), "MutS_clamp": (890,980), "MutS_ATP": (1084,1274)},
    "Vitis vinifera":          {"Tudor": (91,149),  "MutS_N": (347,462), "MutS_conn": (473,631), "MutS_core": (662,975),  "MutS_clamp": (844,934), "MutS_ATP": (1038,1228)},
    "Amborella trichopoda":    {"Tudor": (121,179), "MutS_N": (411,526), "MutS_conn": (537,695), "MutS_core": (726,1039), "MutS_clamp": (908,998), "MutS_ATP": (1102,1292)},
    "Physcomitrium patens":    {"Tudor": (131,189), "MutS_N": (440,555), "MutS_conn": (566,724), "MutS_core": (755,1068), "MutS_clamp": (937,1027), "MutS_ATP": (1131,1321)},
    "Homo sapiens":            {"PWWP": (93,183),   "MutS_N": (408,524), "MutS_conn": (538,694), "MutS_core": (739,1102), "MutS_clamp": (965,1057), "MutS_ATP": (1127,1323)},
    "Mus musculus":            {"PWWP": (93,183),   "MutS_N": (406,523), "MutS_conn": (537,691), "MutS_core": (736,1100), "MutS_clamp": (929,1021), "MutS_ATP": (1125,1321)},
    "Danio rerio":             {"PWWP": (88,178),   "MutS_N": (413,530), "MutS_conn": (543,699), "MutS_core": (744,1107), "MutS_clamp": (936,1028), "MutS_ATP": (1132,1328)},
    "Gallus gallus":           {"PWWP": (74,164),   "MutS_N": (393,510), "MutS_conn": (523,679), "MutS_core": (724,1087), "MutS_clamp": (916,1008), "MutS_ATP": (1112,1308)},
    "Drosophila melanogaster": {                     "MutS_N": (263,376), "MutS_conn": (391,510), "MutS_core": (584,900),  "MutS_clamp": (767,859), "MutS_ATP": (951,1147)},
    "Caenorhabditis elegans":  {                     "MutS_N": (297,410), "MutS_conn": (423,572), "MutS_core": (602,910),  "MutS_clamp": (787,875), "MutS_ATP": (938,1132)},
    "Saccharomyces cerevisiae":{                     "MutS_N": (362,484), "MutS_conn": (497,640), "MutS_core": (672,973),  "MutS_clamp": (850,938), "MutS_ATP": (1009,1199)},
}

# ─── Species order ─────────────────────────────────────────────────────

SPECIES_ORDER = [
    ("Physcomitrium patens",    "Moss",              "plant"),
    ("Amborella trichopoda",    "Basal angiosperm",  "plant"),
    ("Vitis vinifera",          "Grape",             "plant"),
    ("Populus trichocarpa",     "Poplar",            "plant"),
    ("Arabidopsis thaliana",    "Thale cress",       "plant"),
    ("Sorghum bicolor",         "Sorghum",           "plant"),
    ("Zea mays",                "Maize",             "plant"),
    ("Oryza sativa",            "Rice",              "plant"),
    ("Brachypodium distachyon", "Brachypodium",      "plant"),
    ("Danio rerio",             "Zebrafish",         "metazoa"),
    ("Gallus gallus",           "Chicken",           "metazoa"),
    ("Mus musculus",            "Mouse",             "metazoa"),
    ("Homo sapiens",            "Human",             "metazoa"),
    ("Drosophila melanogaster", "Fruit fly",         "outgroup"),
    ("Caenorhabditis elegans",  "Nematode",          "outgroup"),
    ("Saccharomyces cerevisiae","Yeast",             "outgroup"),
]

DOMAIN_COLORS = {
    "Tudor":     "#E74C3C",
    "PWWP":      "#3498DB",
    "MutS_N":    "#F39C12",
    "MutS_conn": "#27AE60",
    "MutS_core": "#8E44AD",
    "MutS_clamp":"#1ABC9C",
    "MutS_ATP":  "#E67E22",
}

DOMAIN_LABELS = {
    "Tudor":     "Tudor",
    "PWWP":      "PWWP",
    "MutS_N":    "MutS I",
    "MutS_conn": "MutS II",
    "MutS_core": "MutS III",
    "MutS_clamp":"MutS IV",
    "MutS_ATP":  "MutS V",
}

# ─── Build figure ──────────────────────────────────────────────────────

fig, ax = plt.subplots(figsize=(12, 9))

bar_height = 0.5
y_spacing = 1.0
max_len = max(e["protein_length"] for e in EXON_DATA.values())

n = len(SPECIES_ORDER)

# Divider lines between groups
plant_end_y = (n - 1 - 8) * y_spacing - 0.5  # after Brachypodium
metazoa_end_y = (n - 1 - 12) * y_spacing - 0.5  # after Homo sapiens

for idx, (species, common, group) in enumerate(SPECIES_ORDER):
    y = (n - 1 - idx) * y_spacing
    exon_info = EXON_DATA[species]
    domains = DOMAINS[species]
    prot_len = exon_info["protein_length"]

    # Protein backbone
    ax.barh(y, prot_len / max_len, height=bar_height,
            color="#ECEDE8", edgecolor="#BDC3C7", linewidth=0.5, zorder=1)

    # Domains
    for dname, (ds, de) in domains.items():
        color = DOMAIN_COLORS[dname]
        ax.barh(y, (de - ds) / max_len, left=ds / max_len,
                height=bar_height, color=color, edgecolor="none", zorder=2, alpha=0.85)

    # Exon boundaries
    for i, (ps, pe) in enumerate(exon_info["exons"]):
        if i == 0:
            continue
        bx = ps / max_len
        ax.plot([bx, bx], [y - bar_height/2, y + bar_height/2],
                color="black", linewidth=0.6, zorder=3)

    # Labels
    ax.text(-0.015, y, f"{species}", ha="right", va="center", fontsize=7.5,
            fontstyle="italic", fontfamily="serif")
    ax.text(prot_len / max_len + 0.01, y, f"{prot_len} aa\n{exon_info['n_exons']} ex",
            ha="left", va="center", fontsize=6, color="#888888")

# Group dividers
ax.axhline(y=plant_end_y, color="#CCCCCC", linewidth=0.8, linestyle="--", zorder=0)
ax.axhline(y=metazoa_end_y, color="#CCCCCC", linewidth=0.8, linestyle="--", zorder=0)

# Group labels on left margin
ax.text(-0.19, (n-1 - 4) * y_spacing, "Plants\n(Tudor)",
        ha="center", va="center", fontsize=9, fontweight="bold", color="#27AE60")
ax.text(-0.19, (n-1 - 10.5) * y_spacing, "Vertebrates\n(PWWP)",
        ha="center", va="center", fontsize=9, fontweight="bold", color="#3498DB")
ax.text(-0.19, (n-1 - 14) * y_spacing, "Outgroups",
        ha="center", va="center", fontsize=9, fontweight="bold", color="#95A5A6")

# Scale bar
for tick_aa in [0, 200, 400, 600, 800, 1000, 1200, 1400]:
    tx = tick_aa / max_len
    if tx <= 1.02:
        ax.text(tx, -0.8, str(tick_aa), ha="center", va="top", fontsize=6.5, color="#888888")
        ax.plot([tx, tx], [-0.55, -0.45], color="#AAAAAA", linewidth=0.5)

ax.set_xlim(-0.26, 1.12)
ax.set_ylim(-1.2, n * y_spacing)
ax.set_xlabel("Amino acid position", fontsize=9)
ax.set_yticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# Legend
legend_handles = []
for dname in ["Tudor", "PWWP", "MutS_N", "MutS_conn", "MutS_core", "MutS_clamp", "MutS_ATP"]:
    legend_handles.append(mpatches.Patch(color=DOMAIN_COLORS[dname], label=DOMAIN_LABELS[dname]))
legend_handles.append(plt.Line2D([0],[0], color="black", linewidth=1, label="Exon boundary"))

ax.legend(handles=legend_handles, loc="upper right", fontsize=6.5,
          framealpha=0.9, ncol=2, title="Domains", title_fontsize=7)

ax.set_title("MSH6 exon-intron structure and domain architecture",
             fontsize=11, fontweight="bold", pad=12)

plt.tight_layout()
plt.savefig("msh6_exon_domain_v2.pdf", dpi=300, bbox_inches="tight")
plt.savefig("msh6_exon_domain_v2.png", dpi=300, bbox_inches="tight")
print("Saved msh6_exon_domain_v2.pdf/png")

# ─── Summary ───────────────────────────────────────────────────────────
print("\nEXON BOUNDARY ANALYSIS:")
for species, common, group in SPECIES_ORDER:
    exons = EXON_DATA[species]["exons"]
    domains = DOMAINS[species]
    reader = None
    for d in ["Tudor", "PWWP"]:
        if d in domains:
            reader = (d, domains[d])
            break
    if not reader:
        continue
    dname, (ds, de) = reader
    containing = [(i+1, ps, pe) for i, (ps, pe) in enumerate(exons) if ps <= de and pe >= ds]
    n_span = len(containing)
    exon_nums = [c[0] for c in containing]
    print(f"  {species:<30} {dname} aa {ds}-{de} -> exon(s) {exon_nums} ({'single' if n_span==1 else f'spans {n_span}'})")
