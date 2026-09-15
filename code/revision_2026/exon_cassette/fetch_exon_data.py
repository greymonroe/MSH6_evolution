#!/usr/bin/env python3
"""
Fetch MSH6 exon-intron structure from NCBI for selected species.
Outputs: exon_data.tsv, domain_data.tsv for R figure construction.

Pipeline: protein accession -> gene ID -> genomic location -> GenBank CDS -> exon boundaries
Domain coords: from all_msh6_interproscan.tsv, keyed to same protein accession.
"""

import csv, json, re, sys, time
from pathlib import Path
from Bio import Entrez, SeqIO
import requests

Entrez.email = 'greymonroe@gmail.com'
WORK = Path(__file__).parent
INTERPRO_FILE = WORK / 'all_msh6_interproscan.tsv'

# Curated species list: (tree_name, protein_accession, reader_type)
# reader_type: 'Tudor', 'PWWP', or 'none'
# Protein accessions chosen from predicted_msh6_blast_annotated.txt
# Prefer RefSeq (NP_/XP_) from GCF_ assemblies
SPECIES = [
    # Plants (Tudor)
    ('Arabidopsis_thaliana',    'NP_192116.1',     'Tudor'),
    ('Amborella_trichopoda',    'XP_006841417.1',  'Tudor'),
    ('Vitis_vinifera',          'XP_002275930.1',  'Tudor'),
    ('Populus_trichocarpa',     'KAI5565098.1',    'Tudor'),   # GCA_ assembly
    ('Oryza_sativa',            'XP_015611811.1',  'Tudor'),
    ('Zea_mays',                'XP_020396805.1',  'Tudor'),
    ('Sorghum_bicolor',         'KXG35601.1',      'Tudor'),   # GCA_ assembly
    ('Physcomitrium_patens',    'PNR32598.1',      'Tudor'),   # GCA_ assembly
    ('Marchantia_paleacea',     'PTQ40051.1',      'Tudor'),   # GCA_ assembly

    # Non-vascular / algae (no reader)
    ('Selaginella_moellendorffii', 'XP_024542030.1', 'none'),
    ('Chlamydomonas_reinhardtii',  'PNW83581.1',     'none'),  # GCA_

    # Sponge (no reader)
    ('Amphimedon_queenslandica', 'XP_019855210.1', 'none'),

    # Cnidaria
    ('Nematostella_vectensis',  'EDO28217.1',      'none'),    # cnidarian without PWWP
    ('Acropora_millepora',      'XP_044184150.1',  'PWWP'),
    ('Hydractinia_symbiolongicarpus', 'XP_057309213.1', 'PWWP'),

    # Protostomes (PWWP)
    ('Daphnia_magna',           'XP_046461715.1',  'PWWP'),
    ('Octopus_sinensis',        'XP_029635395.1',  'PWWP'),
    ('Ixodes_scapularis',       'XP_029829430.3',  'PWWP'),
    ('Capitella_teleta',        'ELU18588.1',      'PWWP'),    # GCA_
    ('Helobdella_robusta',      'ESO12543.1',      'PWWP'),    # GCA_

    # Basal deuterostomes
    ('Strongylocentrotus_purpuratus', 'XP_030829946.1', 'PWWP'),
    ('Branchiostoma_belcheri',  'XP_019632875.1',  'PWWP'),

    # Vertebrates (PWWP)
    ('Petromyzon_marinus',      'XP_032811773.1',  'PWWP'),
    ('Danio_rerio',             'NP_878280.1',     'PWWP'),
    ('Xenopus_tropicalis',      'KAE8600724.1',    'PWWP'),    # GCA_
    ('Anolis_carolinensis',     'XP_008100917.1',  'PWWP'),
    ('Gallus_gallus',           'XP_015133604.2',  'PWWP'),
    ('Mus_musculus',            'NP_034960.1',     'PWWP'),
    ('Homo_sapiens',            'NP_000170.1',     'PWWP'),

    # Outgroups (no reader)
    ('Drosophila_melanogaster', 'NP_648755.1',     'none'),
    ('Caenorhabditis_elegans',  'NP_491163.1',     'none'),
    ('Saccharomyces_cerevisiae','NP_010382.3',     'none'),
]


def get_exon_structure(protein_acc, tree_name):
    """Fetch exon structure from NCBI for a protein accession."""
    print(f'  [{tree_name}] {protein_acc}', flush=True)

    # Step 1: protein -> gene ID
    try:
        handle = Entrez.elink(dbfrom='protein', db='gene', id=protein_acc)
        links = Entrez.read(handle)
        handle.close()
        time.sleep(0.4)
    except Exception as e:
        return None, f'elink failed: {e}'

    if not links[0]['LinkSetDb']:
        return None, 'no gene link'

    gene_id = links[0]['LinkSetDb'][0]['Link'][0]['Id']
    print(f'    gene_id={gene_id}', flush=True)

    # Step 2: gene -> genomic location via NCBI Datasets API
    try:
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}'
        r = requests.get(url, headers={'Accept': 'application/json'}, timeout=30)
        r.raise_for_status()
        time.sleep(0.4)
    except Exception as e:
        return None, f'datasets API failed: {e}'

    report = r.json().get('reports', [{}])[0]
    gene_desc = report.get('gene', {}).get('description', '?')
    gene_symbol = report.get('gene', {}).get('symbol', '?')
    print(f'    gene: {gene_symbol} — {gene_desc}', flush=True)

    annotations = report.get('gene', {}).get('annotations', [])
    if not annotations:
        return None, 'no annotations'

    locs = annotations[0].get('genomic_locations', [])
    if not locs:
        return None, 'no genomic locations'

    # Prefer NC_ (RefSeq chromosome)
    loc = None
    for l in locs:
        gacc = l.get('genomic_accession_version', '')
        if gacc.startswith('NC_'):
            loc = l
            break
    if loc is None:
        loc = locs[0]

    gacc = loc['genomic_accession_version']
    grange = loc['genomic_range']
    start = int(grange['begin'])
    end = int(grange['end'])
    print(f'    genomic: {gacc}:{start}-{end}', flush=True)

    # Step 3: fetch genomic GenBank slice
    try:
        handle = Entrez.efetch(db='nuccore', id=gacc, rettype='gb', retmode='text',
                               seq_start=max(1, start - 2000), seq_stop=end + 2000)
        record = SeqIO.read(handle, 'genbank')
        handle.close()
        time.sleep(0.4)
    except Exception as e:
        return None, f'efetch failed: {e}'

    # Step 4: find CDS matching our exact protein_id
    for feat in record.features:
        if feat.type != 'CDS':
            continue
        pid = feat.qualifiers.get('protein_id', [''])[0]
        if pid != protein_acc:
            continue

        parts = feat.location.parts
        exons = []
        cum_nt = 0
        for i, p in enumerate(parts):
            exon_nt = int(p.end) - int(p.start)
            exon_aa_start = cum_nt // 3 + 1
            cum_nt += exon_nt
            exon_aa_end = cum_nt // 3
            exons.append({
                'exon_num': i + 1,
                'start_aa': exon_aa_start,
                'end_aa': exon_aa_end,
                'length_nt': exon_nt,
                'length_aa': exon_aa_end - exon_aa_start + 1,
            })
        total_aa = cum_nt // 3
        print(f'    matched: {len(exons)} exons, {total_aa} aa', flush=True)

        return {
            'species': tree_name,
            'protein': protein_acc,
            'gene_id': gene_id,
            'gene_symbol': gene_symbol,
            'n_exons': len(exons),
            'total_aa': total_aa,
            'exons': exons,
        }, None

    # CDS not found — list what we did find
    seen = []
    for feat in record.features:
        if feat.type == 'CDS':
            pid = feat.qualifiers.get('protein_id', [''])[0]
            if pid:
                seen.append(pid)
    return None, f'CDS protein_id mismatch. Found: {seen[:8]}'


def get_interproscan_domains(protein_acc):
    """Pull domain annotations from local InterProScan TSV."""
    domain_sigs = {
        'Tudor': {'SM00333', 'cd21828'},  # SMART TUDOR_7, CDD Tudor_Agenet
        'PWWP': {'PF00855', 'SM00293'},
        'MutS_I': {'PF01624'},
        'MutS_II': {'PF05188'},
        'MutS_III': {'PF05192'},
        'MutS_IV': {'PF05190'},
        'MutS_V': {'PF00488'},
    }

    domains = []
    with open(INTERPRO_FILE) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            if cols[0] != protein_acc:
                continue
            sig_acc = cols[4]
            start = int(cols[6])
            end = int(cols[7])
            for domain_name, sigs in domain_sigs.items():
                if sig_acc in sigs:
                    domains.append({
                        'domain': domain_name,
                        'sig_acc': sig_acc,
                        'start_aa': start,
                        'end_aa': end,
                    })
    return domains


def main():
    exon_rows = []
    domain_rows = []
    summary_rows = []
    failures = []

    print(f'Fetching exon structure for {len(SPECIES)} species...\n', flush=True)

    for tree_name, protein_acc, reader_type in SPECIES:
        result, error = get_exon_structure(protein_acc, tree_name)

        if error:
            print(f'    *** FAILED: {error}', flush=True)
            failures.append((tree_name, protein_acc, error))
            continue

        # Exon data
        for ex in result['exons']:
            exon_rows.append({
                'species': tree_name,
                'protein': protein_acc,
                'reader_type': reader_type,
                'exon_num': ex['exon_num'],
                'start_aa': ex['start_aa'],
                'end_aa': ex['end_aa'],
                'length_aa': ex['length_aa'],
                'length_nt': ex['length_nt'],
                'total_aa': result['total_aa'],
                'n_exons': result['n_exons'],
            })

        # Domain data from InterProScan
        domains = get_interproscan_domains(protein_acc)
        for dom in domains:
            domain_rows.append({
                'species': tree_name,
                'protein': protein_acc,
                'reader_type': reader_type,
                'domain': dom['domain'],
                'sig_acc': dom['sig_acc'],
                'start_aa': dom['start_aa'],
                'end_aa': dom['end_aa'],
                'total_aa': result['total_aa'],
            })

        summary_rows.append({
            'species': tree_name,
            'protein': protein_acc,
            'reader_type': reader_type,
            'n_exons': result['n_exons'],
            'total_aa': result['total_aa'],
            'gene_symbol': result['gene_symbol'],
            'n_domains': len(domains),
        })

        print(f'    domains: {[d["domain"] for d in domains]}', flush=True)
        print(flush=True)

    # Write exon TSV
    exon_path = WORK / 'exon_data.tsv'
    with open(exon_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['species', 'protein', 'reader_type', 'exon_num',
                           'start_aa', 'end_aa', 'length_aa', 'length_nt', 'total_aa', 'n_exons'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(exon_rows)
    print(f'\nWrote {len(exon_rows)} exon rows to {exon_path}')

    # Write domain TSV
    domain_path = WORK / 'domain_data.tsv'
    with open(domain_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['species', 'protein', 'reader_type', 'domain',
                           'sig_acc', 'start_aa', 'end_aa', 'total_aa'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(domain_rows)
    print(f'Wrote {len(domain_rows)} domain rows to {domain_path}')

    # Write summary
    summary_path = WORK / 'species_summary.tsv'
    with open(summary_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['species', 'protein', 'reader_type', 'n_exons',
                           'total_aa', 'gene_symbol', 'n_domains'],
                           delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)
    print(f'Wrote {len(summary_rows)} species to {summary_path}')

    if failures:
        print(f'\n*** {len(failures)} failures:')
        for name, acc, err in failures:
            print(f'  {name} ({acc}): {err}')


if __name__ == '__main__':
    main()
