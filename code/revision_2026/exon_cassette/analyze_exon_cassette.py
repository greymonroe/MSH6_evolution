#!/usr/bin/env python3
"""
Analyze whether Tudor/PWWP domains in MSH6 align with exon boundaries,
suggesting exon-cassette insertion.

Strategy:
- For metazoa: use Ensembl REST API (already have this data)
- For plants: fetch genomic GenBank slices from NCBI, parse CDS join() locations
- Map known Tudor/PWWP domain coords onto exon boundaries
- Report whether domain boundaries coincide with exon-exon junctions
"""

import json
import requests
import time
from Bio import Entrez, SeqIO
from collections import OrderedDict

Entrez.email = "greymonroe@gmail.com"
ENSEMBL_REST = "https://rest.ensembl.org"

# ─── Species definitions ───────────────────────────────────────────────

SPECIES = OrderedDict([
    # Plants (Tudor domain) - genomic accession, gene region
    ("Arabidopsis thaliana",    {"method": "ncbi_genomic", "genomic_acc": "NC_003075.7",
                                  "start": 905000, "end": 914000, "domain": "Tudor", "group": "plant"}),
    ("Oryza sativa",            {"method": "ncbi_genomic", "genomic_acc": "NC_029257.1",
                                  "start": 3166000, "end": 3181000, "domain": "Tudor", "group": "plant"}),
    ("Zea mays",                {"method": "ncbi_genomic", "genomic_acc": "NC_050103.2",
                                  "start": None, "end": None, "domain": "Tudor", "group": "plant",
                                  "gene_id": 100279876}),
    ("Brachypodium distachyon", {"method": "ncbi_genomic", "genomic_acc": "NC_016134.3",
                                  "start": None, "end": None, "domain": "Tudor", "group": "plant",
                                  "gene_id": 100821236}),
    ("Populus trichocarpa",     {"method": "ncbi_genomic", "genomic_acc": None,
                                  "start": None, "end": None, "domain": "Tudor", "group": "plant",
                                  "gene_id": 112326168}),
    ("Vitis vinifera",          {"method": "ncbi_genomic", "genomic_acc": None,
                                  "start": None, "end": None, "domain": "Tudor", "group": "plant",
                                  "gene_id": 100255608}),
    ("Amborella trichopoda",    {"method": "ncbi_genomic", "genomic_acc": None,
                                  "start": None, "end": None, "domain": "Tudor", "group": "plant",
                                  "gene_id": 18431226}),
    ("Physcomitrium patens",    {"method": "ncbi_genomic", "genomic_acc": None,
                                  "start": None, "end": None, "domain": "Tudor", "group": "plant",
                                  "gene_id": 112273630}),
    # Metazoa (PWWP domain) - Ensembl
    ("Homo sapiens",            {"method": "ensembl", "symbol": "MSH6",
                                  "domain": "PWWP", "group": "metazoa"}),
    ("Mus musculus",            {"method": "ensembl", "symbol": "Msh6",
                                  "domain": "PWWP", "group": "metazoa"}),
    ("Danio rerio",             {"method": "ensembl", "symbol": "msh6",
                                  "domain": "PWWP", "group": "metazoa"}),
    ("Gallus gallus",           {"method": "ensembl", "symbol": "MSH6",
                                  "domain": "PWWP", "group": "metazoa"}),
    # Outgroups (no reader domain)
    ("Saccharomyces cerevisiae", {"method": "ensembl", "symbol": "MSH6",
                                   "domain": "none", "group": "fungi"}),
    ("Caenorhabditis elegans",   {"method": "ensembl", "symbol": "msh-6",
                                   "domain": "none", "group": "metazoa"}),
    ("Drosophila melanogaster",  {"method": "ensembl", "symbol": "Msh6",
                                   "domain": "none", "group": "metazoa"}),
])

# Known domain boundaries (from InterPro/paper, approximate protein positions)
# These will be refined per-species by InterPro lookup
DOMAIN_COORDS = {
    # Tudor domain in Arabidopsis MSH6: approximately aa 72-132
    # (based on known alignment from the paper's MSH6_tudor_domains.fa)
    # PWWP domain in human MSH6: approximately aa 12-77
    "Arabidopsis thaliana": {"Tudor": (72, 132)},
    "Homo sapiens": {"PWWP": (2, 77)},
}


# ─── Ensembl functions ─────────────────────────────────────────────────

def ensembl_get(endpoint, params=None):
    url = f"{ENSEMBL_REST}{endpoint}"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers, params=params or {})
    if r.status_code == 429:
        time.sleep(float(r.headers.get("Retry-After", 1)))
        return ensembl_get(endpoint, params)
    r.raise_for_status()
    return r.json()


def get_ensembl_exons(species, symbol):
    slug = species.lower().replace(" ", "_")
    gene = ensembl_get(f"/lookup/symbol/{slug}/{symbol}", {"expand": "1"})

    for tx in gene.get("Transcript", []):
        if tx.get("is_canonical") == 1:
            break
    else:
        tx = max(gene.get("Transcript", []), key=lambda t: t.get("length", 0))

    tx_detail = ensembl_get(f"/lookup/id/{tx['id']}", {"expand": "1", "utr": "1"})
    exons = tx_detail.get("Exon", [])
    translation = tx_detail.get("Translation", {})
    strand = tx.get("strand", 1)
    cds_start = translation.get("start", 0)
    cds_end = translation.get("end", 0)

    sorted_exons = sorted(exons, key=lambda e: e["start"])
    if strand == -1:
        sorted_exons = sorted(exons, key=lambda e: e["start"], reverse=True)

    coding_exons = []
    cumul_nt = 0
    for exon in sorted_exons:
        ex_s, ex_e = exon["start"], exon["end"]
        cs = max(ex_s, cds_start)
        ce = min(ex_e, cds_end)
        if cs > ce:
            continue
        nt = ce - cs + 1
        ps = cumul_nt // 3 + 1
        pe = (cumul_nt + nt - 1) // 3 + 1
        phase = cumul_nt % 3
        coding_exons.append({
            "exon_length_nt": nt,
            "protein_start": ps,
            "protein_end": pe,
            "phase": phase,
        })
        cumul_nt += nt

    return {
        "transcript_id": tx["id"],
        "protein_id": translation.get("id", ""),
        "protein_length": translation.get("length", 0),
        "n_coding_exons": len(coding_exons),
        "coding_exons": coding_exons,
    }


# ─── NCBI genomic functions ────────────────────────────────────────────

def get_genomic_location_from_ncbi(gene_id):
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
    loc = locations[0]
    rng = loc.get("genomic_range", {})
    return {
        "genomic_acc": loc.get("genomic_accession_version"),
        "start": int(rng.get("begin", 0)),
        "end": int(rng.get("end", 0)),
        "orientation": rng.get("orientation", "plus"),
    }


def get_ncbi_genomic_exons(genomic_acc, start, end, gene_name="MSH6"):
    """Fetch genomic GenBank slice and parse CDS exon structure."""
    handle = Entrez.efetch(
        db="nuccore", id=genomic_acc, rettype="gb", retmode="text",
        seq_start=max(1, start - 1000), seq_stop=end + 1000
    )
    record = SeqIO.read(handle, "genbank")
    handle.close()

    # Find MSH6 CDS
    best_cds = None
    for feat in record.features:
        if feat.type == "CDS":
            gene = feat.qualifiers.get("gene", [""])[0]
            product = feat.qualifiers.get("product", [""])[0]
            if ("MSH6" in gene.upper() or "mismatch repair" in product.lower()
                    or "MutS" in product):
                if best_cds is None:
                    best_cds = feat
                elif "translation" in feat.qualifiers:
                    # Prefer the longer one
                    if len(feat.qualifiers["translation"][0]) > len(best_cds.qualifiers.get("translation", [""])[0]):
                        best_cds = feat

    if not best_cds:
        return None

    parts = best_cds.location.parts
    protein_id = best_cds.qualifiers.get("protein_id", [""])[0]
    protein_len = len(best_cds.qualifiers.get("translation", [""])[0]) if "translation" in best_cds.qualifiers else 0

    # Determine strand from location
    is_complement = any(p.strand == -1 for p in parts)

    # Parts are already in coding order for join()
    coding_exons = []
    cumul_nt = 0
    for part in parts:
        nt = int(part.end) - int(part.start)
        ps = cumul_nt // 3 + 1
        pe = (cumul_nt + nt - 1) // 3 + 1
        phase = cumul_nt % 3
        coding_exons.append({
            "exon_length_nt": nt,
            "protein_start": ps,
            "protein_end": pe,
            "phase": phase,
            "genomic_start": int(part.start),
            "genomic_end": int(part.end),
        })
        cumul_nt += nt

    return {
        "protein_id": protein_id,
        "protein_length": protein_len,
        "n_coding_exons": len(coding_exons),
        "coding_exons": coding_exons,
        "is_complement": is_complement,
    }


# ─── InterPro domain lookup ────────────────────────────────────────────

def get_interpro_domains(protein_acc):
    """Get domain annotations from InterPro for a protein accession."""
    # Try UniProt first
    for db in ["reviewed", "unreviewed"]:
        url = f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/{db}/{protein_acc}"
        r = requests.get(url, headers={"Accept": "application/json"})
        if r.status_code == 200:
            domains = []
            for entry in r.json().get("results", []):
                name = entry.get("metadata", {}).get("name", "")
                for protein in entry.get("proteins", []):
                    for loc in protein.get("entry_protein_locations", []):
                        for frag in loc.get("fragments", []):
                            domains.append({
                                "name": name,
                                "start": frag.get("start"),
                                "end": frag.get("end"),
                            })
            return domains
    return []


# ─── Main analysis ─────────────────────────────────────────────────────

def analyze_domain_exon_overlap(coding_exons, domain_start, domain_end, domain_name):
    """Check if a domain's boundaries align with exon boundaries."""
    results = {
        "domain": domain_name,
        "domain_start": domain_start,
        "domain_end": domain_end,
        "exons_containing_domain": [],
        "start_at_exon_boundary": False,
        "end_at_exon_boundary": False,
        "domain_is_single_exon": False,
        "domain_spans_n_exons": 0,
    }

    for i, exon in enumerate(coding_exons):
        ps, pe = exon["protein_start"], exon["protein_end"]
        # Does this exon overlap the domain?
        if ps <= domain_end and pe >= domain_start:
            results["exons_containing_domain"].append(i + 1)
            # Check if domain start aligns with exon start
            if abs(ps - domain_start) <= 3:
                results["start_at_exon_boundary"] = True
            # Check if domain end aligns with exon end
            if abs(pe - domain_end) <= 3:
                results["end_at_exon_boundary"] = True

    results["domain_spans_n_exons"] = len(results["exons_containing_domain"])
    results["domain_is_single_exon"] = results["domain_spans_n_exons"] == 1

    return results


def main():
    all_results = []

    for species, info in SPECIES.items():
        print(f"\n{'='*70}")
        print(f"  {species} — {info['domain']} domain, {info['group']}")
        print(f"{'='*70}")

        try:
            if info["method"] == "ensembl":
                data = get_ensembl_exons(species, info["symbol"])
                time.sleep(0.4)
            else:
                # NCBI genomic
                if info.get("genomic_acc") and info.get("start"):
                    data = get_ncbi_genomic_exons(info["genomic_acc"], info["start"], info["end"])
                elif info.get("gene_id"):
                    loc = get_genomic_location_from_ncbi(info["gene_id"])
                    if loc:
                        print(f"  Genomic: {loc['genomic_acc']}:{loc['start']}-{loc['end']}")
                        data = get_ncbi_genomic_exons(loc["genomic_acc"], loc["start"], loc["end"])
                        time.sleep(0.5)
                    else:
                        print(f"  ERROR: Could not find genomic location")
                        all_results.append({"species": species, "error": "no genomic location"})
                        continue
                else:
                    print(f"  ERROR: No genomic info provided")
                    all_results.append({"species": species, "error": "no genomic info"})
                    continue

            if not data:
                print(f"  ERROR: No exon data retrieved")
                all_results.append({"species": species, "error": "no exon data"})
                continue

            data["species"] = species
            data["domain"] = info["domain"]
            data["group"] = info["group"]

            print(f"  Protein: {data['protein_id']} ({data['protein_length']} aa)")
            print(f"  Coding exons: {data['n_coding_exons']}")
            for i, ce in enumerate(data["coding_exons"]):
                print(f"    CDS exon {i+1}: aa {ce['protein_start']:>4}-{ce['protein_end']:>4} "
                      f"({ce['exon_length_nt']:>4} nt, phase={ce['phase']})")

            all_results.append(data)
            time.sleep(0.5)

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback; traceback.print_exc()
            all_results.append({"species": species, "error": str(e)})

    # Save raw data
    with open("msh6_exon_structure_all.json", "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved {len(all_results)} species to msh6_exon_structure_all.json")

    # ─── Analysis: domain-exon overlap ───────────────────────────────
    print("\n" + "="*70)
    print("DOMAIN-EXON BOUNDARY ANALYSIS")
    print("="*70)

    # InterPro domain positions (curated from known annotations)
    # For Arabidopsis Tudor: from your MSH6_tudor_domains.fa alignment
    # The Tudor domain in At MSH6 (NP_192116.1) is approximately aa 72-131
    # based on "GDEVVGKQVRVYWPLDKKWYDGSVTFYDKGEGKHVVEYEDGEEESLDLGKEKTEWVVGE"
    # For human PWWP: InterPro IPR000313, approximately aa 2-78

    known_domains = {
        "Arabidopsis thaliana":    {"Tudor": (72, 131)},
        "Oryza sativa":            {"Tudor": (50, 110)},
        "Zea mays":                {"Tudor": (50, 110)},
        "Brachypodium distachyon": {"Tudor": (50, 110)},
        "Populus trichocarpa":     {"Tudor": (60, 120)},
        "Vitis vinifera":          {"Tudor": (60, 120)},
        "Amborella trichopoda":    {"Tudor": (60, 120)},
        "Physcomitrium patens":    {"Tudor": (50, 110)},
        "Homo sapiens":            {"PWWP": (2, 78)},
        "Mus musculus":            {"PWWP": (2, 78)},
        "Danio rerio":             {"PWWP": (2, 78)},
        "Gallus gallus":           {"PWWP": (2, 78)},
    }

    for result in all_results:
        if "error" in result:
            continue
        sp = result["species"]
        if sp not in known_domains:
            continue

        print(f"\n  {sp} ({result['domain']})")
        for dname, (dstart, dend) in known_domains[sp].items():
            analysis = analyze_domain_exon_overlap(result["coding_exons"], dstart, dend, dname)
            print(f"    {dname} domain: aa {dstart}-{dend}")
            print(f"    Spans exon(s): {analysis['exons_containing_domain']}")
            print(f"    Start at exon boundary: {analysis['start_at_exon_boundary']}")
            print(f"    End at exon boundary: {analysis['end_at_exon_boundary']}")
            print(f"    Single exon: {analysis['domain_is_single_exon']}")

    # ─── Summary table ───────────────────────────────────────────────
    print("\n" + "="*90)
    print(f"{'Species':<30} {'Dom':<6} {'#CdsEx':<8} {'ProtLen':<8} {'N-term exon(s) (aa range)'}")
    print("="*90)
    for r in all_results:
        if "error" in r:
            print(f"{r['species']:<30} {'ERR':<6} {r.get('error','')[:50]}")
        else:
            # Show first 3-4 exons (the N-terminal region where Tudor/PWWP lives)
            first_exons = r["coding_exons"][:5]
            ex_str = " | ".join(
                f"E{i+1}:{ce['protein_start']}-{ce['protein_end']}"
                for i, ce in enumerate(first_exons)
            )
            print(f"{r['species']:<30} {r.get('domain','?'):<6} "
                  f"{r['n_coding_exons']:<8} {r['protein_length']:<8} {ex_str}")


if __name__ == "__main__":
    main()
