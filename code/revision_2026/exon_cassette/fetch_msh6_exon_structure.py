#!/usr/bin/env python3
"""
Fetch MSH6 exon-intron structure from Ensembl REST API and map
Tudor/PWWP domain boundaries onto exon coordinates.

Uses:
- Ensembl REST (ensembl.org) for metazoa/fungi
- Ensembl Plants REST (rest.ensembl.org with plants division) for plants
- InterPro API for domain coordinates on the protein
"""

import json
import requests
import time
import sys

ENSEMBL_REST = "https://rest.ensembl.org"


def ensembl_get(endpoint, params=None):
    """Generic Ensembl REST API GET."""
    url = f"{ENSEMBL_REST}{endpoint}"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers, params=params or {})
    if r.status_code == 429:
        wait = float(r.headers.get("Retry-After", 1))
        time.sleep(wait)
        return ensembl_get(endpoint, params)
    r.raise_for_status()
    return r.json()


def get_gene_by_symbol(species, symbol):
    """Look up gene by symbol."""
    species_slug = species.lower().replace(" ", "_")
    try:
        data = ensembl_get(f"/lookup/symbol/{species_slug}/{symbol}",
                           {"expand": "1"})
        return data
    except Exception as e:
        print(f"    Lookup failed for {symbol} in {species}: {e}")
        return None


def get_transcript_exons(transcript_id):
    """Get exon details for a transcript."""
    data = ensembl_get(f"/lookup/id/{transcript_id}",
                       {"expand": "1", "utr": "1"})
    return data


def get_interpro_domains(uniprot_acc):
    """Get domain annotations from InterPro for a UniProt accession."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{uniprot_acc}"
    r = requests.get(url, headers={"Accept": "application/json"})
    if r.status_code != 200:
        return []
    data = r.json()
    domains = []
    for entry in data.get("results", []):
        name = entry.get("metadata", {}).get("name", "")
        accession = entry.get("metadata", {}).get("accession", "")
        entry_type = entry.get("metadata", {}).get("type", "")
        for protein in entry.get("proteins", []):
            for loc in protein.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    domains.append({
                        "name": name,
                        "accession": accession,
                        "type": entry_type,
                        "start": frag.get("start"),
                        "end": frag.get("end"),
                    })
    return domains


def get_interpro_domains_by_protein(protein_acc):
    """Get domain annotations from InterPro by protein accession (UniProt or other)."""
    url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/reviewed/{protein_acc}"
    r = requests.get(url, headers={"Accept": "application/json"})
    if r.status_code != 200:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/unreviewed/{protein_acc}"
        r = requests.get(url, headers={"Accept": "application/json"})
    if r.status_code != 200:
        return []
    return parse_interpro_response(r.json())


def parse_interpro_response(data):
    domains = []
    for entry in data.get("results", []):
        name = entry.get("metadata", {}).get("name", "")
        accession = entry.get("metadata", {}).get("accession", "")
        for protein in entry.get("proteins", []):
            for loc in protein.get("entry_protein_locations", []):
                for frag in loc.get("fragments", []):
                    domains.append({
                        "name": name,
                        "accession": accession,
                        "start": frag.get("start"),
                        "end": frag.get("end"),
                    })
    return domains


SPECIES = [
    # Plants with Tudor domain in MSH6
    {"name": "Arabidopsis thaliana",    "symbol": "AT4G02070", "domain": "Tudor",
     "uniprot": "F4JLG3", "group": "plant"},
    {"name": "Oryza sativa",            "symbol": "Os02g0142200", "domain": "Tudor",
     "uniprot": None, "group": "plant"},
    {"name": "Zea mays",                "symbol": "Zm00001eb092860", "domain": "Tudor",
     "uniprot": None, "group": "plant"},
    {"name": "Brachypodium distachyon", "symbol": "BRADI_3g04670", "domain": "Tudor",
     "uniprot": None, "group": "plant"},
    {"name": "Populus trichocarpa",     "symbol": "Potri.014G070500", "domain": "Tudor",
     "uniprot": None, "group": "plant"},
    {"name": "Vitis vinifera",          "symbol": "GSVIVG01010131001", "domain": "Tudor",
     "uniprot": None, "group": "plant"},
    # Metazoa with PWWP domain in MSH6
    {"name": "Homo sapiens",            "symbol": "MSH6",  "domain": "PWWP",
     "uniprot": "P52701", "group": "metazoa"},
    {"name": "Mus musculus",            "symbol": "Msh6",  "domain": "PWWP",
     "uniprot": "P54276", "group": "metazoa"},
    {"name": "Danio rerio",             "symbol": "msh6",  "domain": "PWWP",
     "uniprot": None, "group": "metazoa"},
    {"name": "Gallus gallus",           "symbol": "MSH6",  "domain": "PWWP",
     "uniprot": None, "group": "metazoa"},
    # Outgroups (no reader domain)
    {"name": "Saccharomyces cerevisiae", "symbol": "MSH6", "domain": "none",
     "uniprot": "P54820", "group": "fungi"},
    {"name": "Caenorhabditis elegans",   "symbol": "msh-6", "domain": "none",
     "uniprot": None, "group": "metazoa"},
    {"name": "Drosophila melanogaster",  "symbol": "Msh6",  "domain": "none",
     "uniprot": None, "group": "metazoa"},
]


def compute_exon_protein_coords(exons, cds_start_genomic, cds_end_genomic, strand):
    """
    Map each exon to protein coordinates based on CDS boundaries.
    exons: list of dicts with 'start', 'end' (1-based genomic coords)
    cds_start_genomic, cds_end_genomic: CDS boundaries (1-based)
    strand: 1 or -1
    """
    if strand == 1:
        sorted_exons = sorted(exons, key=lambda e: e["start"])
    else:
        sorted_exons = sorted(exons, key=lambda e: e["start"], reverse=True)

    results = []
    cumul_coding_nt = 0

    for exon in sorted_exons:
        ex_start = exon["start"]
        ex_end = exon["end"]

        # Clip to CDS
        if strand == 1:
            coding_start = max(ex_start, cds_start_genomic)
            coding_end = min(ex_end, cds_end_genomic)
        else:
            coding_start = max(ex_start, cds_start_genomic)
            coding_end = min(ex_end, cds_end_genomic)

        if coding_start > coding_end:
            results.append({
                "exon_start": ex_start, "exon_end": ex_end,
                "exon_length_nt": ex_end - ex_start + 1,
                "coding_nt": 0,
                "protein_start": None, "protein_end": None,
                "phase": None, "is_utr": True,
            })
            continue

        coding_nt = coding_end - coding_start + 1
        prot_start = cumul_coding_nt // 3 + 1
        prot_end = (cumul_coding_nt + coding_nt - 1) // 3 + 1
        phase = cumul_coding_nt % 3

        results.append({
            "exon_start": ex_start, "exon_end": ex_end,
            "exon_length_nt": ex_end - ex_start + 1,
            "coding_nt": coding_nt,
            "protein_start": prot_start,
            "protein_end": prot_end,
            "phase": phase,
            "is_utr": False,
        })
        cumul_coding_nt += coding_nt

    return results


def process_species(sp):
    """Process one species: get gene model, exon structure, domain locations."""
    name = sp["name"]
    symbol = sp["symbol"]
    print(f"\n{'='*60}")
    print(f"  {name} — symbol={symbol}, expected domain={sp['domain']}")
    print(f"{'='*60}")

    # Step 1: Look up gene
    gene = get_gene_by_symbol(name, symbol)
    if not gene:
        return {"species": name, "error": "gene lookup failed"}

    gene_id = gene.get("id", "")
    print(f"  Gene: {gene_id} ({gene.get('display_name', '')})")
    print(f"  Biotype: {gene.get('biotype', '')}")

    # Step 2: Get canonical transcript
    transcripts = gene.get("Transcript", [])
    if not transcripts:
        return {"species": name, "error": "no transcripts"}

    # Prefer canonical, else longest
    canonical = None
    for tx in transcripts:
        if tx.get("is_canonical", 0) == 1:
            canonical = tx
            break
    if not canonical:
        canonical = max(transcripts, key=lambda t: t.get("length", 0))

    tx_id = canonical["id"]
    print(f"  Canonical transcript: {tx_id}")
    time.sleep(0.35)

    # Step 3: Get full exon info
    tx_detail = get_transcript_exons(tx_id)
    exons = tx_detail.get("Exon", [])
    print(f"  Exons: {len(exons)}")

    # Translation (CDS) info
    translation = tx_detail.get("Translation")
    if not translation:
        return {"species": name, "error": "no translation/CDS", "n_exons": len(exons)}

    protein_id = translation.get("id", "")
    protein_length = translation.get("length", 0)
    cds_start = translation.get("start", 0)
    cds_end = translation.get("end", 0)
    print(f"  Protein: {protein_id} ({protein_length} aa)")
    print(f"  CDS genomic: {cds_start}-{cds_end}")

    strand = canonical.get("strand", 1)

    # Step 4: Compute protein coordinates per exon
    exon_data = [{"start": e["start"], "end": e["end"]} for e in exons]
    exon_map = compute_exon_protein_coords(exon_data, cds_start, cds_end, strand)

    coding_exons = [e for e in exon_map if not e["is_utr"]]
    print(f"  Coding exons: {len(coding_exons)}")
    for i, ce in enumerate(coding_exons):
        print(f"    CDS exon {i+1}: aa {ce['protein_start']}-{ce['protein_end']} "
              f"({ce['coding_nt']} nt, phase={ce['phase']})")

    time.sleep(0.35)

    return {
        "species": name,
        "domain": sp["domain"],
        "group": sp["group"],
        "gene_id": gene_id,
        "display_name": gene.get("display_name", ""),
        "transcript_id": tx_id,
        "protein_id": protein_id,
        "protein_length": protein_length,
        "strand": strand,
        "n_exons_total": len(exons),
        "n_coding_exons": len(coding_exons),
        "coding_exons": coding_exons,
    }


def main():
    all_results = []

    for sp in SPECIES:
        try:
            result = process_species(sp)
            all_results.append(result)
        except Exception as e:
            print(f"  ERROR: {e}")
            all_results.append({"species": sp["name"], "error": str(e)})
        time.sleep(0.5)

    # Save
    outfile = "msh6_exon_structure.json"
    with open(outfile, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved to {outfile}")

    # Summary table
    print(f"\n{'='*90}")
    print(f"{'Species':<30} {'Domain':<8} {'Protein':<20} {'CodExons':<10} {'ProtLen':<10}")
    print(f"{'='*90}")
    for r in all_results:
        if "error" in r:
            print(f"{r['species']:<30} {'ERR':<8} {r.get('error','')[:40]}")
        else:
            print(f"{r['species']:<30} {r['domain']:<8} {r['protein_id']:<20} "
                  f"{r['n_coding_exons']:<10} {r['protein_length']:<10}")


if __name__ == "__main__":
    main()
