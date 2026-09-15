#!/usr/bin/env python3
"""
Fetch MSH6 exon structure for plant species from NCBI RefSeq mRNA GenBank records,
and combine with the Ensembl-fetched metazoa data.

For plants: fetch RefSeq mRNA -> get linked genomic record -> parse exon structure
Alternative: use the NCBI Datasets API gene download endpoint.
"""

from Bio import Entrez, SeqIO
import json
import time

Entrez.email = "greymonroe@gmail.com"

# RefSeq mRNA accessions for plant MSH6 (one per species, preferred longest/canonical)
PLANT_MRNAS = {
    "Arabidopsis thaliana":    {"acc": "NM_116441.4",    "domain": "Tudor", "gene_id": 828147},
    "Oryza sativa":            {"acc": "XM_015756325.3", "domain": "Tudor", "gene_id": 4347016},
    "Zea mays":                {"acc": "NM_001152829.1", "domain": "Tudor", "gene_id": 100279876},
    "Brachypodium distachyon": {"acc": "XM_003576466.4", "domain": "Tudor", "gene_id": 100821236},
    "Populus trichocarpa":     {"acc": "XM_052447133.1", "domain": "Tudor", "gene_id": 112326168},
    "Vitis vinifera":          {"acc": "XM_002275894.4", "domain": "Tudor", "gene_id": 100255608},
    "Amborella trichopoda":    {"acc": "XM_006841354.3", "domain": "Tudor", "gene_id": 18431226},
    "Physcomitrium patens":    {"acc": "XM_024502604.2", "domain": "Tudor", "gene_id": 112273630},
}


def fetch_mrna_genbank(accession):
    """Fetch GenBank record for a RefSeq mRNA."""
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def get_genomic_location(gene_id):
    """Get genomic accession and coordinates from NCBI Datasets API."""
    import requests
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gene_id}"
    r = requests.get(url, headers={"Accept": "application/json"})
    if r.status_code != 200:
        return None
    data = r.json()
    reports = data.get("reports", [])
    if not reports:
        return None
    gene = reports[0].get("gene", {})
    annotations = gene.get("annotations", [])
    if not annotations:
        return None
    locations = annotations[0].get("genomic_locations", [])
    if not locations:
        return None
    return {
        "genomic_acc": locations[0].get("genomic_accession_version"),
        "seq_name": locations[0].get("sequence_name"),
        "start": int(locations[0].get("genomic_range", {}).get("begin", 0)),
        "end": int(locations[0].get("genomic_range", {}).get("end", 0)),
        "orientation": locations[0].get("genomic_range", {}).get("orientation", "plus"),
    }


def fetch_genomic_genbank(genomic_acc, start, end):
    """Fetch a slice of a genomic GenBank record."""
    handle = Entrez.efetch(
        db="nuccore", id=genomic_acc, rettype="gb", retmode="text",
        seq_start=start, seq_stop=end
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record


def parse_mrna_features(record):
    """Parse CDS and other features from mRNA GenBank record."""
    info = {
        "accession": record.id,
        "length": len(record.seq),
        "description": record.description,
        "cds": None,
        "protein_id": None,
        "protein_length": None,
        "coded_by": None,
    }

    for feat in record.features:
        if feat.type == "CDS":
            info["cds"] = {
                "start": int(feat.location.start),
                "end": int(feat.location.end),
            }
            if "protein_id" in feat.qualifiers:
                info["protein_id"] = feat.qualifiers["protein_id"][0]
            if "translation" in feat.qualifiers:
                info["protein_length"] = len(feat.qualifiers["translation"][0])
                info["protein_seq"] = feat.qualifiers["translation"][0]
            if "coded_by" in feat.qualifiers:
                info["coded_by"] = feat.qualifiers["coded_by"][0]

    return info


def parse_genomic_gene_features(record, gene_id=None):
    """Parse exon/CDS/mRNA features from a genomic GenBank record for the MSH6 gene."""
    genes = []
    mrnas = []
    cds_features = []
    exon_features = []

    for feat in record.features:
        if feat.type == "gene":
            genes.append(feat)
        elif feat.type == "mRNA":
            mrnas.append(feat)
        elif feat.type == "CDS":
            cds_features.append(feat)
        elif feat.type == "exon":
            exon_features.append(feat)

    # Find the MSH6 gene
    msh6_gene = None
    for g in genes:
        note = " ".join(g.qualifiers.get("note", []))
        gene_name = g.qualifiers.get("gene", [""])[0]
        db_xref = g.qualifiers.get("db_xref", [])
        if "MSH6" in gene_name.upper() or "mismatch repair" in note.lower():
            msh6_gene = g
            break
        if gene_id:
            for xref in db_xref:
                if str(gene_id) in xref:
                    msh6_gene = g
                    break

    return {
        "n_genes": len(genes),
        "n_mrnas": len(mrnas),
        "n_cds": len(cds_features),
        "n_exons": len(exon_features),
        "msh6_gene": str(msh6_gene.location) if msh6_gene else None,
        "mrnas": [{"location": str(m.location), "product": m.qualifiers.get("product", [""])[0]}
                  for m in mrnas[:3]],
        "cds": [{"location": str(c.location), "protein_id": c.qualifiers.get("protein_id", [""])[0]}
                for c in cds_features[:3]],
    }


def get_exon_structure_from_coded_by(coded_by_str):
    """
    Parse the 'coded_by' qualifier from a CDS feature on an mRNA record.
    Format: "join(NC_xxx:1234..5678,NC_xxx:7890..9012,...)"
    Returns list of exon coords on the genomic sequence.
    """
    if not coded_by_str:
        return []

    import re
    # Handle complement() wrapper
    is_complement = "complement" in coded_by_str
    inner = coded_by_str
    if "complement" in inner:
        inner = re.sub(r"complement\((.*)\)", r"\1", inner)

    # Handle join()
    if "join" in inner:
        inner = re.sub(r"join\((.*)\)", r"\1", inner)

    parts = inner.split(",")
    exons = []
    for part in parts:
        part = part.strip()
        # Parse "ACC:start..end" or "start..end"
        match = re.match(r"(?:([^:]+):)?<?(\d+)\.\.>?(\d+)", part)
        if match:
            acc = match.group(1) or ""
            start = int(match.group(2))
            end = int(match.group(3))
            exons.append({"acc": acc, "start": start, "end": end})

    return exons, is_complement


def compute_protein_coords_from_cds_exons(exons):
    """Given CDS exon ranges (genomic), compute protein coordinate per exon."""
    results = []
    cumul_nt = 0
    for exon in exons:
        nt_len = exon["end"] - exon["start"] + 1
        prot_start = cumul_nt // 3 + 1
        prot_end = (cumul_nt + nt_len - 1) // 3 + 1
        phase = cumul_nt % 3
        results.append({
            "genomic_start": exon["start"],
            "genomic_end": exon["end"],
            "exon_length_nt": nt_len,
            "protein_start": prot_start,
            "protein_end": prot_end,
            "phase": phase,
        })
        cumul_nt += nt_len
    return results


def main():
    all_results = []

    # Load existing metazoa results
    try:
        with open("msh6_exon_structure.json") as f:
            metazoa_results = json.load(f)
        print(f"Loaded {len(metazoa_results)} existing results from Ensembl")
        all_results = metazoa_results
    except FileNotFoundError:
        print("No existing results found")

    print("\n" + "="*70)
    print("Fetching plant MSH6 exon structures from NCBI")
    print("="*70)

    plant_results = []

    for species, info in PLANT_MRNAS.items():
        print(f"\n{'—'*60}")
        print(f"  {species} — mRNA: {info['acc']}")
        print(f"{'—'*60}")

        try:
            # Fetch mRNA GenBank record
            record = fetch_mrna_genbank(info["acc"])
            mrna_info = parse_mrna_features(record)
            print(f"  mRNA length: {mrna_info['length']} nt")
            print(f"  Protein: {mrna_info['protein_id']} ({mrna_info['protein_length']} aa)")
            print(f"  CDS on mRNA: {mrna_info['cds']['start']}-{mrna_info['cds']['end']}")

            coded_by = mrna_info.get("coded_by")
            if coded_by:
                print(f"  coded_by: {coded_by[:100]}...")
                exons, is_complement = get_exon_structure_from_coded_by(coded_by)
                if is_complement:
                    exons = list(reversed(exons))
                print(f"  CDS exons from coded_by: {len(exons)}")

                exon_prot = compute_protein_coords_from_cds_exons(exons)
                for i, ep in enumerate(exon_prot):
                    print(f"    CDS exon {i+1}: aa {ep['protein_start']}-{ep['protein_end']} "
                          f"({ep['exon_length_nt']} nt, phase={ep['phase']})")

                result = {
                    "species": species,
                    "domain": info["domain"],
                    "group": "plant",
                    "gene_id": str(info["gene_id"]),
                    "transcript_id": info["acc"],
                    "protein_id": mrna_info["protein_id"],
                    "protein_length": mrna_info["protein_length"],
                    "n_coding_exons": len(exons),
                    "coding_exons": exon_prot,
                    "coded_by": coded_by,
                }
                plant_results.append(result)
            else:
                # No coded_by — try to get exon structure from genomic record
                print("  No coded_by qualifier — trying genomic approach...")
                genomic = get_genomic_location(info["gene_id"])
                if genomic:
                    print(f"  Genomic: {genomic['genomic_acc']}:{genomic['start']}-{genomic['end']}")
                    # For now just note it
                    result = {
                        "species": species,
                        "domain": info["domain"],
                        "group": "plant",
                        "gene_id": str(info["gene_id"]),
                        "transcript_id": info["acc"],
                        "protein_id": mrna_info["protein_id"],
                        "protein_length": mrna_info["protein_length"],
                        "genomic_location": genomic,
                        "note": "no coded_by, need genomic GFF parse",
                    }
                    plant_results.append(result)
                else:
                    plant_results.append({
                        "species": species, "error": "no exon structure available",
                        "protein_id": mrna_info["protein_id"],
                        "protein_length": mrna_info["protein_length"],
                    })

        except Exception as e:
            print(f"  ERROR: {e}")
            plant_results.append({"species": species, "error": str(e)})

        time.sleep(0.5)

    # Combine
    all_results.extend(plant_results)

    # Save
    outfile = "msh6_exon_structure_all.json"
    with open(outfile, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved {len(all_results)} species to {outfile}")

    # Summary table
    print(f"\n{'='*90}")
    print(f"{'Species':<30} {'Domain':<8} {'Protein':<20} {'CodExons':<10} {'ProtLen':<10}")
    print(f"{'='*90}")
    for r in all_results:
        if "error" in r:
            print(f"{r['species']:<30} {'ERR':<8} {r.get('error','')[:40]}")
        else:
            print(f"{r['species']:<30} {r.get('domain','?'):<8} "
                  f"{r.get('protein_id',''):<20} "
                  f"{r.get('n_coding_exons','?'):<10} "
                  f"{r.get('protein_length','?'):<10}")


if __name__ == "__main__":
    main()
