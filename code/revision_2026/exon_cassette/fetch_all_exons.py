#!/usr/bin/env python3
"""
Fetch MSH6 exon-intron structure for ALL single-copy MSH6 organisms.
Two pipelines:
  GenBank proteins -> efetch coded_by (batch, 1 call per batch)
  RefSeq proteins  -> elink gene -> Datasets API -> efetch genomic CDS (3 calls each)

Outputs: exon_data_all.tsv, domain_data_all.tsv, species_summary_all.tsv
"""

import csv, json, re, sys, time, traceback
from pathlib import Path
from Bio import Entrez, SeqIO
from io import StringIO
import requests

Entrez.email = 'greymonroe@gmail.com'
WORK = Path(__file__).parent
INTERPRO_FILE = WORK / 'all_msh6_interproscan.tsv'
BLAST_FILE = WORK / 'predicted_msh6_blast_annotated.txt'
ORG_FILE = WORK / 'msh6_organism_table.txt'

DOMAIN_SIGS = {
    'Tudor': {'SM00333', 'cd20404', 'cd04508', 'cd20401'},
    'PWWP': {'PF00855', 'SM00293', 'cd05837'},
    'MutS_I': {'PF01624'}, 'MutS_II': {'PF05188'}, 'MutS_III': {'PF05192'},
    'MutS_IV': {'PF05190'}, 'MutS_V': {'PF00488'},
}


def load_organism_table():
    orgs = {}
    with open(ORG_FILE) as f:
        f.readline()
        for line in f:
            m = re.match(r'"([^"]+)"(.*)', line.strip())
            if not m:
                continue
            org = m.group(1)
            fields = m.group(2).strip().split()
            orgs[org] = {
                'muts_count': int(fields[0]),
                'tudor_call': fields[5] if len(fields) >= 6 else 'Absent',
                'pwwp_call': fields[6] if len(fields) >= 7 else 'Absent',
            }
    return orgs


def load_blast_table():
    rows = []
    with open(BLAST_FILE) as f:
        f.readline()
        for line in f:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            gc = parts[0]
            protein = parts[1]
            m = re.search(r'"([^"]+)"', line)
            if not m:
                continue
            org = m.group(1)
            after = line[m.end():].strip().split()
            muts = after[0] if after else 'FALSE'
            tudor_interpro = after[-3] if len(after) >= 4 else 'FALSE'
            pwwp_interpro = after[-1] if len(after) >= 1 else 'FALSE'
            if muts == 'TRUE':
                rows.append({
                    'gc': gc, 'protein': protein, 'org': org,
                    'tudor': tudor_interpro == 'TRUE',
                    'pwwp': pwwp_interpro == 'TRUE',
                })
    return rows


def select_representatives(org_data, blast_rows):
    """Pick one protein per organism. Prefer RefSeq, then reader-domain-bearing."""
    by_org = {}
    for row in blast_rows:
        org = row['org']
        if org not in by_org:
            by_org[org] = []
        by_org[org].append(row)

    selected = {}
    for org, info in org_data.items():
        if org not in by_org:
            continue
        candidates = by_org[org]

        reader_type = 'none'
        if info['tudor_call'] == 'Present':
            reader_type = 'Tudor'
        elif info['pwwp_call'] == 'Present':
            reader_type = 'PWWP'

        best = None
        for c in candidates:
            is_refseq = c['protein'].startswith(('XP_', 'NP_'))
            has_reader = c['tudor'] or c['pwwp']
            score = (is_refseq * 2) + (has_reader * 1)
            if best is None or score > best[1]:
                best = (c, score)

        if best:
            c = best[0]
            tree_name = org.replace(' ', '_')
            selected[tree_name] = {
                'protein': c['protein'],
                'gc': c['gc'],
                'reader_type': reader_type,
                'is_refseq': c['protein'].startswith(('XP_', 'NP_')),
            }
    return selected


def parse_coded_by(coded_by_str):
    parts = re.findall(r'(\d+)\.\.(\d+)', coded_by_str)
    if not parts:
        return None
    exons = []
    cum_nt = 0
    for i, (s, e) in enumerate(parts):
        exon_nt = int(e) - int(s) + 1
        aa_s = cum_nt // 3 + 1
        cum_nt += exon_nt
        aa_e = cum_nt // 3
        exons.append({
            'exon_num': i + 1, 'start_aa': aa_s, 'end_aa': aa_e,
            'length_aa': aa_e - aa_s + 1, 'length_nt': exon_nt,
        })
    return exons, cum_nt // 3


def batch_fetch_coded_by(protein_accs, batch_size=150):
    """Batch fetch protein records and extract coded_by exon structures."""
    results = {}
    for i in range(0, len(protein_accs), batch_size):
        batch = protein_accs[i:i + batch_size]
        batch_num = i // batch_size + 1
        total_batches = (len(protein_accs) + batch_size - 1) // batch_size
        print(f'  Batch {batch_num}/{total_batches} ({len(batch)} proteins)...', flush=True)

        try:
            handle = Entrez.efetch(db='protein', id=','.join(batch),
                                   rettype='gp', retmode='text')
            text = handle.read()
            handle.close()
            time.sleep(0.5)
        except Exception as e:
            print(f'    Batch failed: {e}', flush=True)
            time.sleep(2)
            continue

        for record in SeqIO.parse(StringIO(text), 'genbank'):
            acc = record.id
            if '.' not in acc:
                acc = record.name
            for feat in record.features:
                if feat.type == 'CDS':
                    coded_by = feat.qualifiers.get('coded_by', [None])[0]
                    if coded_by and 'join' in coded_by:
                        result = parse_coded_by(coded_by)
                        if result:
                            exons, total_aa = result
                            results[acc] = {
                                'n_exons': len(exons),
                                'total_aa': total_aa,
                                'exons': exons,
                            }
                            break
                    elif coded_by:
                        result = parse_coded_by(coded_by)
                        if result:
                            exons, total_aa = result
                            results[acc] = {
                                'n_exons': len(exons),
                                'total_aa': total_aa,
                                'exons': exons,
                            }
                            break

    return results


def fetch_refseq_exons(protein_acc):
    """Fetch exon structure for RefSeq protein via gene pipeline."""
    try:
        handle = Entrez.elink(dbfrom='protein', db='gene', id=protein_acc)
        links = Entrez.read(handle)
        handle.close()
        time.sleep(0.35)
    except Exception:
        return None

    if not links[0]['LinkSetDb']:
        return None
    gene_id = links[0]['LinkSetDb'][0]['Link'][0]['Id']

    try:
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}'
        r = requests.get(url, headers={'Accept': 'application/json'}, timeout=30)
        r.raise_for_status()
        time.sleep(0.35)
    except Exception:
        return None

    try:
        report = r.json()['reports'][0]
        locs = report['gene']['annotations'][0]['genomic_locations']
    except (KeyError, IndexError):
        return None

    loc = None
    for l in locs:
        if l.get('genomic_accession_version', '').startswith('NC_'):
            loc = l
            break
    if not loc:
        loc = locs[0]

    gacc = loc['genomic_accession_version']
    start = int(loc['genomic_range']['begin'])
    end = int(loc['genomic_range']['end'])

    try:
        handle = Entrez.efetch(db='nuccore', id=gacc, rettype='gb', retmode='text',
                               seq_start=max(1, start - 2000), seq_stop=end + 2000)
        record = SeqIO.read(handle, 'genbank')
        handle.close()
        time.sleep(0.35)
    except Exception:
        return None

    for feat in record.features:
        if feat.type != 'CDS':
            continue
        pid = feat.qualifiers.get('protein_id', [''])[0]
        if pid != protein_acc:
            continue
        parts = feat.location.parts
        exons = []
        cum = 0
        for i, p in enumerate(parts):
            nt = int(p.end) - int(p.start)
            s = cum // 3 + 1
            cum += nt
            e = cum // 3
            exons.append({
                'exon_num': i + 1, 'start_aa': s, 'end_aa': e,
                'length_aa': e - s + 1, 'length_nt': nt,
            })
        return {'n_exons': len(exons), 'total_aa': cum // 3, 'exons': exons}
    return None


def load_interproscan_index():
    """Build index: protein_acc -> list of domain hits."""
    print('Loading InterProScan index...', flush=True)
    index = {}
    with open(INTERPRO_FILE) as f:
        for line in f:
            cols = line.strip().split('\t')
            if len(cols) < 9:
                continue
            acc = cols[0]
            sig_acc = cols[4]
            for dname, sigs in DOMAIN_SIGS.items():
                if sig_acc in sigs:
                    if acc not in index:
                        index[acc] = []
                    index[acc].append({
                        'domain': dname, 'sig_acc': sig_acc,
                        'start_aa': int(cols[6]), 'end_aa': int(cols[7]),
                    })
    print(f'  Indexed {len(index)} proteins with domain hits', flush=True)
    return index


def main():
    print('=== Loading data ===', flush=True)
    org_data = load_organism_table()
    blast_rows = load_blast_table()
    selected = select_representatives(org_data, blast_rows)
    interpro = load_interproscan_index()

    print(f'Selected {len(selected)} representative proteins', flush=True)

    gb_proteins = {k: v for k, v in selected.items() if not v['is_refseq']}
    rs_proteins = {k: v for k, v in selected.items() if v['is_refseq']}
    print(f'  GenBank: {len(gb_proteins)}')
    print(f'  RefSeq:  {len(rs_proteins)}')

    # Phase 1: batch coded_by for GenBank proteins
    print(f'\n=== Phase 1: GenBank coded_by ({len(gb_proteins)} proteins) ===', flush=True)
    gb_accs = [v['protein'] for v in gb_proteins.values()]
    coded_by_results = batch_fetch_coded_by(gb_accs)
    print(f'  Got exon data for {len(coded_by_results)} / {len(gb_accs)} GenBank proteins', flush=True)

    # Phase 2: gene pipeline for RefSeq proteins
    print(f'\n=== Phase 2: RefSeq gene pipeline ({len(rs_proteins)} proteins) ===', flush=True)
    rs_results = {}
    for i, (sp, info) in enumerate(rs_proteins.items()):
        if (i + 1) % 50 == 0:
            print(f'  {i+1}/{len(rs_proteins)}...', flush=True)
        result = fetch_refseq_exons(info['protein'])
        if result:
            rs_results[info['protein']] = result
    print(f'  Got exon data for {len(rs_results)} / {len(rs_proteins)} RefSeq proteins', flush=True)

    # Combine and write
    all_results = {**coded_by_results, **rs_results}
    print(f'\n=== Writing output ({len(all_results)} proteins total) ===', flush=True)

    exon_rows = []
    domain_rows = []
    summary_rows = []

    for sp, info in sorted(selected.items()):
        protein = info['protein']
        if protein not in all_results:
            continue
        result = all_results[protein]
        reader_type = info['reader_type']

        for ex in result['exons']:
            exon_rows.append({
                'species': sp, 'protein': protein, 'reader_type': reader_type,
                'exon_num': ex['exon_num'], 'start_aa': ex['start_aa'],
                'end_aa': ex['end_aa'], 'length_aa': ex['length_aa'],
                'length_nt': ex['length_nt'], 'total_aa': result['total_aa'],
                'n_exons': result['n_exons'],
            })

        domains = interpro.get(protein, [])
        for d in domains:
            domain_rows.append({
                'species': sp, 'protein': protein, 'reader_type': reader_type,
                'domain': d['domain'], 'sig_acc': d['sig_acc'],
                'start_aa': d['start_aa'], 'end_aa': d['end_aa'],
                'total_aa': result['total_aa'],
            })

        summary_rows.append({
            'species': sp, 'protein': protein, 'reader_type': reader_type,
            'gc': info['gc'], 'n_exons': result['n_exons'],
            'total_aa': result['total_aa'], 'n_domains': len(domains),
        })

    exon_path = WORK / 'exon_data_all.tsv'
    with open(exon_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'species', 'protein', 'reader_type', 'exon_num', 'start_aa',
            'end_aa', 'length_aa', 'length_nt', 'total_aa', 'n_exons',
        ], delimiter='\t')
        w.writeheader()
        w.writerows(exon_rows)

    domain_path = WORK / 'domain_data_all.tsv'
    with open(domain_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'species', 'protein', 'reader_type', 'domain', 'sig_acc',
            'start_aa', 'end_aa', 'total_aa',
        ], delimiter='\t')
        w.writeheader()
        w.writerows(domain_rows)

    summary_path = WORK / 'species_summary_all.tsv'
    with open(summary_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'species', 'protein', 'reader_type', 'gc', 'n_exons',
            'total_aa', 'n_domains',
        ], delimiter='\t')
        w.writeheader()
        w.writerows(summary_rows)

    print(f'\nWrote {len(exon_rows)} exon rows to {exon_path}')
    print(f'Wrote {len(domain_rows)} domain rows to {domain_path}')
    print(f'Wrote {len(summary_rows)} species to {summary_path}')

    # Stats
    n_tudor = sum(1 for r in summary_rows if r['reader_type'] == 'Tudor')
    n_pwwp = sum(1 for r in summary_rows if r['reader_type'] == 'PWWP')
    n_none = sum(1 for r in summary_rows if r['reader_type'] == 'none')
    print(f'\n  Tudor: {n_tudor}, PWWP: {n_pwwp}, none: {n_none}')


if __name__ == '__main__':
    main()
