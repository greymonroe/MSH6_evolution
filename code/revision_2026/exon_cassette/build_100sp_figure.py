#!/usr/bin/env python3
"""
Build combined dataset for ~100-species figure.
Merge curated 41-species data with comprehensive dataset, prune TimeTree.
"""

import csv, re
from pathlib import Path
from Bio import Phylo

WORK = Path(__file__).parent

def load_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))

# Load both datasets
curated_exons = load_tsv(WORK / 'exon_data.tsv')
curated_doms  = load_tsv(WORK / 'domain_data.tsv')
curated_summ  = load_tsv(WORK / 'species_summary.tsv')
all_exons     = load_tsv(WORK / 'exon_data_all.tsv')
all_doms      = load_tsv(WORK / 'domain_data_all.tsv')
all_summ      = load_tsv(WORK / 'species_summary_all.tsv')

curated_species = {r['species'] for r in curated_summ}
all_species_set = {r['species'] for r in all_summ}

# TimeTree
tree = Phylo.read('/Users/greymonroe/repos/tol_reader_repair/files/time_tree_may4.nwk', 'newick')
tree_tips = {t.name for t in tree.get_terminals()}
tree_genus = {}
for t in tree_tips:
    p = t.split('_')
    if len(p) >= 2:
        tree_genus[p[0] + '_' + p[1]] = t

def find_tree_name(sp):
    if sp in tree_tips:
        return sp
    p = sp.split('_')
    if len(p) >= 2:
        key = p[0] + '_' + p[1]
        if key in tree_genus:
            return tree_genus[key]
    return None

# Species to include from comprehensive dataset
comp_targets = [
    'Capitella_teleta', 'Mizuhopecten_yessoensis',
    'Centruroides_sculpturatus', 'Limulus_polyphemus',
    'Callorhinchus_milii', 'Latimeria_chalumnae',
    'Aedes_aegypti', 'Tenebrio_molitor', 'Hypsibius_exemplaris',
    'Neurospora_crassa', 'Aspergillus_niger', 'Cryptococcus_neoformans',
    'Candida_albicans', 'Schizosaccharomyces_pombe',
    'Plasmodium_falciparum_3D7', 'Toxoplasma_gondii',
    'Tetrahymena_thermophila_SB210', 'Paramecium_tetraurelia',
    'Chondrus_crispus', 'Porphyra_umbilicalis',
    'Glycine_max', 'Solanum_lycopersicum', 'Helianthus_annuus',
    'Coffea_arabica', 'Citrus_sinensis', 'Theobroma_cacao',
    'Asparagus_officinalis', 'Nymphaea_colorata',
    'Beta_vulgaris_subsp._vulgaris', 'Cannabis_sativa',
    'Chlorella_vulgaris', 'Coccomyxa_subellipsoidea_C-169',
    'Chara_braunii', 'Klebsormidium_nitens',
    'Micromonas_pusilla_CCMP1545', 'Diphasiastrum_complanatum',
    'Zostera_marina', 'Daphnia_pulex', 'Phytophthora_infestans',
    'Crassostrea_angulata', 'Bombyx_mandarina', 'Apis_mellifera',
    'Trichoplax_sp._H2',
]

def find_comp(target):
    if target in all_species_set:
        return target
    for sp in all_species_set:
        if target.lower() in sp.lower():
            return sp
    return None

# Build final species list
# Start with all 41 curated species
final_species = {}  # data_name -> tree_name
for sp in curated_species:
    tn = find_tree_name(sp)
    if tn:
        final_species[sp] = tn

# Add comprehensive species
for target in comp_targets:
    match = find_comp(target)
    if match:
        tn = find_tree_name(match)
        if tn and match not in final_species:
            final_species[match] = tn

print(f'Final species count: {len(final_species)}')

# Build rename map: tree_name -> data_name
tree_to_data = {v: k for k, v in final_species.items()}

# Rename tree tips to data names
renames = {}
for tip in tree.get_terminals():
    if tip.name in tree_to_data:
        if tip.name != tree_to_data[tip.name]:
            renames[tip.name] = tree_to_data[tip.name]

for tip in tree.get_terminals():
    if tip.name in renames:
        tip.name = renames[tip.name]

# Prune tree
keep = set(final_species.keys())
to_remove = [t.name for t in tree.get_terminals() if t.name not in keep]
for name in to_remove:
    tree.prune(name)

# Remove any remaining non-target tips
extra = {t.name for t in tree.get_terminals()} - keep
for name in extra:
    tree.prune(name)

final_tips = sorted(t.name for t in tree.get_terminals())
print(f'Tree tips: {len(final_tips)}')

# Merge data: curated takes priority, then comprehensive
exon_rows = []
domain_rows = []
summary_rows = []

for sp in final_tips:
    # Check curated first
    if sp in curated_species:
        exon_rows.extend(r for r in curated_exons if r['species'] == sp)
        domain_rows.extend(r for r in curated_doms if r['species'] == sp)
        summary_rows.extend(r for r in curated_summ if r['species'] == sp)
    elif sp in all_species_set:
        exon_rows.extend(r for r in all_exons if r['species'] == sp)
        domain_rows.extend(r for r in all_doms if r['species'] == sp)
        summary_rows.extend(r for r in all_summ if r['species'] == sp)
    else:
        print(f'  WARNING: {sp} not in any dataset')

# Write outputs
def write_tsv(path, rows, fieldnames):
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)

write_tsv(WORK / 'exon_data_100.tsv', exon_rows,
          ['species', 'protein', 'reader_type', 'exon_num', 'start_aa', 'end_aa',
           'length_aa', 'length_nt', 'total_aa', 'n_exons'])

write_tsv(WORK / 'domain_data_100.tsv', domain_rows,
          ['species', 'protein', 'reader_type', 'domain', 'sig_acc',
           'start_aa', 'end_aa', 'total_aa'])

write_tsv(WORK / 'species_summary_100.tsv', summary_rows,
          ['species', 'protein', 'reader_type', 'n_exons', 'total_aa'])

Phylo.write(tree, str(WORK / 'species_tree_100.nwk'), 'newick')

print(f'\nWrote {len(exon_rows)} exon rows, {len(domain_rows)} domain rows')
print(f'{len(summary_rows)} species in summary')
print(f'Tree: species_tree_100.nwk')

# Stats
from collections import Counter
types = Counter(r.get('reader_type', 'none') for r in summary_rows)
print(f'Tudor: {types.get("Tudor",0)}, PWWP: {types.get("PWWP",0)}, none: {types.get("none",0)}')
