#!/usr/bin/env python3
"""
Comprehensive PWWP donor screen.
Strategy:
1. Download all PWWP proteins from InterPro (PF00855) — gets UniProt accessions + domain coords
2. Use UniProt ID mapping to get EMBL/GenBank accessions (these HAVE coded_by)
3. Batch-fetch GenPept records for coded_by exon structure + PWWP sequences
4. Compute exon break positions

Also fetch RefSeq PWWP proteins from NCBI directly (broader coverage).
"""

import csv, json, re, sys, time
from pathlib import Path
from io import StringIO
from collections import defaultdict
from Bio import Entrez, SeqIO
import requests

Entrez.email = 'greymonroe@gmail.com'
WORK = Path(__file__).parent


def fetch_interpro_list():
    """Download all PF00855 proteins from InterPro with domain coords."""
    cache = WORK / 'interpro_pwwp_proteins.json'
    if cache.exists():
        print('  Loading InterPro cache...', flush=True)
        with open(cache) as f:
            return json.load(f)

    proteins = []
    url = 'https://www.ebi.ac.uk/interpro/api/protein/UniProt/entry/pfam/PF00855/?page_size=200&format=json'
    page = 0

    while url:
        page += 1
        if page % 25 == 1:
            print(f'  InterPro page {page} ({len(proteins)} proteins)...', flush=True)
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            data = r.json()
        except Exception as e:
            print(f'    Page {page} failed: {e}', flush=True)
            time.sleep(2)
            continue

        for entry in data.get('results', []):
            meta = entry.get('metadata', {})
            acc = meta.get('accession', '')
            org = meta.get('source_organism', {})
            for loc_group in entry.get('entries', []):
                for loc in loc_group.get('entry_protein_locations', []):
                    for frag in loc.get('fragments', []):
                        proteins.append({
                            'uniprot': acc,
                            'organism': org.get('scientificName', '?'),
                            'tax_id': org.get('taxId', 0),
                            'prot_len': meta.get('length', 0),
                            'pwwp_start': frag.get('start', 0),
                            'pwwp_end': frag.get('end', 0),
                        })
        url = data.get('next')
        time.sleep(0.2)

    with open(cache, 'w') as f:
        json.dump(proteins, f)
    return proteins


def uniprot_to_embl(uniprot_accs, batch_size=500):
    """Map UniProt accessions to EMBL/GenBank protein accessions via UniProt API."""
    mapping = {}
    for i in range(0, len(uniprot_accs), batch_size):
        batch = uniprot_accs[i:i + batch_size]
        batch_num = i // batch_size + 1
        total = (len(uniprot_accs) + batch_size - 1) // batch_size
        if batch_num % 10 == 1:
            print(f'  UniProt mapping batch {batch_num}/{total}...', flush=True)

        try:
            # Submit job
            r = requests.post('https://rest.uniprot.org/idmapping/run',
                              data={'from': 'UniProtKB_AC-ID', 'to': 'EMBL-GenBank-DDBJ_CDS',
                                    'ids': ','.join(batch)}, timeout=30)
            r.raise_for_status()
            job_id = r.json()['jobId']
            time.sleep(1)

            # Poll for results
            for _ in range(30):
                r = requests.get(f'https://rest.uniprot.org/idmapping/status/{job_id}', timeout=30)
                status = r.json()
                if 'results' in status or 'failedIds' in status:
                    break
                if status.get('jobStatus') == 'FINISHED':
                    break
                time.sleep(1)

            # Get results
            r = requests.get(f'https://rest.uniprot.org/idmapping/results/{job_id}?size=500', timeout=30)
            r.raise_for_status()
            results = r.json().get('results', [])
            for res in results:
                uni = res.get('from', '')
                embl = res.get('to', '')
                if uni and embl:
                    if uni not in mapping:
                        mapping[uni] = embl
            time.sleep(0.3)

        except Exception as e:
            print(f'    Mapping batch failed: {e}', flush=True)
            time.sleep(2)

    return mapping


def parse_coded_by(coded_by_str):
    parts = re.findall(r'(\d+)\.\.(\d+)', coded_by_str)
    if not parts:
        return None
    exons = []
    cum_nt = 0
    for i, (s, e) in enumerate(parts):
        nt = int(e) - int(s) + 1
        if nt < 0:
            return None
        aa_s = cum_nt // 3 + 1
        cum_nt += nt
        aa_e = cum_nt // 3
        exons.append({'num': i + 1, 'start': aa_s, 'end': aa_e, 'len': aa_e - aa_s + 1})
    total_aa = cum_nt // 3
    if total_aa < 100:
        return None
    return exons, total_aa


def batch_fetch_ncbi(accessions, batch_size=200):
    """Batch fetch protein records: coded_by + PWWP annotation + sequence."""
    results = {}
    n_batches = (len(accessions) + batch_size - 1) // batch_size

    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        bn = i // batch_size + 1
        if bn % 10 == 1:
            print(f'  NCBI batch {bn}/{n_batches} ({len(results)} results)...', flush=True)

        try:
            handle = Entrez.efetch(db='protein', id=','.join(batch),
                                   rettype='gp', retmode='text')
            text = handle.read()
            handle.close()
            time.sleep(0.5)
        except Exception as e:
            print(f'    Batch {bn} failed: {e}', flush=True)
            time.sleep(2)
            continue

        try:
            for record in SeqIO.parse(StringIO(text), 'genbank'):
                acc = record.id
                seq = str(record.seq)

                coded_by = None
                gene_name = ''
                product = ''
                for feat in record.features:
                    if feat.type == 'CDS':
                        coded_by = feat.qualifiers.get('coded_by', [None])[0]
                        gene_name = feat.qualifiers.get('gene', [''])[0]
                        product = feat.qualifiers.get('product', [''])[0]
                        break

                exons = None
                if coded_by and 'join' in coded_by:
                    result = parse_coded_by(coded_by)
                    if result:
                        exons, _ = result

                pwwp_s = pwwp_e = None
                for feat in record.features:
                    if feat.type == 'Region':
                        name = feat.qualifiers.get('region_name', [''])[0]
                        if 'PWWP' in name.upper():
                            s = int(feat.location.start) + 1
                            e = int(feat.location.end)
                            if pwwp_s is None or s < pwwp_s:
                                pwwp_s = s
                            if pwwp_e is None or e > pwwp_e:
                                pwwp_e = e

                pwwp_seq = ''
                if pwwp_s and pwwp_e and len(seq) >= pwwp_e:
                    pwwp_seq = seq[pwwp_s - 1:pwwp_e]

                results[acc] = {
                    'exons': exons,
                    'prot_len': len(seq),
                    'pwwp_s': pwwp_s, 'pwwp_e': pwwp_e,
                    'pwwp_seq': pwwp_seq,
                    'gene': gene_name, 'product': product,
                }
        except Exception as e:
            print(f'    Parse error: {e}', flush=True)

    return results


def main():
    # ── Phase 1: Get InterPro protein list ──
    print('=== Phase 1: InterPro PWWP protein list ===', flush=True)
    interpro = fetch_interpro_list()
    print(f'  Total entries: {len(interpro)}', flush=True)

    # Deduplicate by UniProt accession
    by_uni = {}
    for p in interpro:
        acc = p['uniprot']
        if acc not in by_uni:
            by_uni[acc] = p
        else:
            if p['pwwp_start'] < by_uni[acc]['pwwp_start']:
                by_uni[acc]['pwwp_start'] = p['pwwp_start']
            if p['pwwp_end'] > by_uni[acc]['pwwp_end']:
                by_uni[acc]['pwwp_end'] = p['pwwp_end']

    print(f'  Unique UniProt accessions: {len(by_uni)}', flush=True)

    # ── Phase 2: Map UniProt -> EMBL/GenBank accessions ──
    print('\n=== Phase 2: UniProt -> EMBL mapping ===', flush=True)
    cache_map = WORK / 'uniprot_to_embl_map.json'
    if cache_map.exists():
        print('  Loading mapping cache...', flush=True)
        with open(cache_map) as f:
            uni_to_embl = json.load(f)
    else:
        uni_accs = list(by_uni.keys())
        uni_to_embl = uniprot_to_embl(uni_accs)
        with open(cache_map, 'w') as f:
            json.dump(uni_to_embl, f)
    print(f'  Mapped {len(uni_to_embl)} / {len(by_uni)} UniProt -> EMBL', flush=True)

    # ── Phase 3: Batch fetch from NCBI ──
    print('\n=== Phase 3: NCBI batch fetch ===', flush=True)
    embl_accs = list(set(uni_to_embl.values()))
    print(f'  Unique EMBL accessions to fetch: {len(embl_accs)}', flush=True)

    cache_ncbi = WORK / 'ncbi_pwwp_results.json'
    if cache_ncbi.exists():
        print('  Loading NCBI cache...', flush=True)
        with open(cache_ncbi) as f:
            ncbi_results = json.load(f)
    else:
        ncbi_results = batch_fetch_ncbi(embl_accs)
        # Serialize (exons are dicts, should be fine)
        with open(cache_ncbi, 'w') as f:
            json.dump(ncbi_results, f)
    print(f'  Got NCBI data for {len(ncbi_results)} proteins', flush=True)

    # ── Phase 4: Compute breaks ──
    print('\n=== Phase 4: Computing PWWP exon breaks ===', flush=True)

    output_rows = []
    pwwp_seqs = []

    for uni_acc, info in by_uni.items():
        embl_acc = uni_to_embl.get(uni_acc)
        if not embl_acc:
            continue
        ncbi = ncbi_results.get(embl_acc)
        if not ncbi:
            continue

        # Use NCBI PWWP annotation if available, else InterPro
        pwwp_s = ncbi.get('pwwp_s') or info['pwwp_start']
        pwwp_e = ncbi.get('pwwp_e') or info['pwwp_end']

        exons = ncbi.get('exons')
        has_exons = exons is not None

        break_rel = None
        n_exons_pwwp = 0
        if exons and pwwp_s and pwwp_e:
            overlapping = [e for e in exons if e['end'] >= pwwp_s and e['start'] <= pwwp_e]
            n_exons_pwwp = len(overlapping)
            if n_exons_pwwp >= 2:
                break_rel = overlapping[0]['end'] - pwwp_s

        pwwp_seq = ncbi.get('pwwp_seq', '')

        output_rows.append({
            'uniprot': uni_acc,
            'embl': embl_acc,
            'organism': info['organism'],
            'tax_id': info['tax_id'],
            'prot_len': ncbi.get('prot_len', info['prot_len']),
            'gene': ncbi.get('gene', ''),
            'product': ncbi.get('product', ''),
            'pwwp_start': pwwp_s,
            'pwwp_end': pwwp_e,
            'pwwp_len': pwwp_e - pwwp_s + 1 if pwwp_s and pwwp_e else 0,
            'has_genomic_exons': has_exons,
            'n_total_exons': len(exons) if exons else 1,
            'n_exons_in_pwwp': n_exons_pwwp,
            'break_from_pwwp_start': break_rel if break_rel is not None else '',
        })

        if pwwp_seq:
            brk = f'break+{break_rel}' if break_rel is not None else 'no_break'
            pwwp_seqs.append(f'>{embl_acc}|{uni_acc} {info["organism"]} {ncbi.get("gene","")} {brk}\n{pwwp_seq}')

    # Write outputs
    out_tsv = WORK / 'pwwp_all_proteins_breaks.tsv'
    with open(out_tsv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=[
            'uniprot', 'embl', 'organism', 'tax_id', 'prot_len', 'gene', 'product',
            'pwwp_start', 'pwwp_end', 'pwwp_len',
            'has_genomic_exons', 'n_total_exons', 'n_exons_in_pwwp',
            'break_from_pwwp_start',
        ], delimiter='\t')
        w.writeheader()
        w.writerows(output_rows)

    out_fasta = WORK / 'pwwp_domain_sequences.fasta'
    with open(out_fasta, 'w') as f:
        f.write('\n'.join(pwwp_seqs))

    print(f'\nWrote {len(output_rows)} rows to {out_tsv}')
    print(f'Wrote {len(pwwp_seqs)} PWWP sequences to {out_fasta}')

    # Summary
    with_exons = [r for r in output_rows if r['has_genomic_exons']]
    with_break = [r for r in output_rows if r['break_from_pwwp_start'] != '']
    print(f'\nWith genomic exon data: {len(with_exons)}')
    print(f'With PWWP exon break: {len(with_break)}')

    if with_break:
        near_62 = [r for r in with_break if abs(int(r['break_from_pwwp_start']) - 62) < 8]
        print(f'\n*** Proteins with break near +62 aa (±8): {len(near_62)} ***')
        for r in sorted(near_62, key=lambda x: abs(int(x['break_from_pwwp_start']) - 62)):
            print(f'  {r["embl"]:20s} {r["gene"]:12s} {r["organism"]:40s} +{r["break_from_pwwp_start"]}aa  {r["product"][:40]}')

        # Distribution of break positions
        from collections import Counter
        breaks = Counter(int(r['break_from_pwwp_start']) for r in with_break)
        print(f'\nBreak position distribution (top 20):')
        for pos, count in breaks.most_common(20):
            bar = '#' * min(count, 50)
            print(f'  +{pos:3d} aa: {count:5d} {bar}')


if __name__ == '__main__':
    main()
