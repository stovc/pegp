#!/usr/bin/env python3
"""
Build a flat list of genome file URLs from GTDB accessions + NCBI assembly_summary files.

Inputs (defaults expect files in genome_info/):
- bac120_metadata_filtered.tsv, ar53_metadata_filtered.tsv (must contain column: accession)
- bacteria_refseq_assembly_summary.txt
- bacteria_genbank_assembly_summary.txt
- archaea_refseq_assembly_summary.txt
- archaea_genbank_assembly_summary.txt

Rules:
- accession looks like RS_GCF_xxx or GB_GCA_xxx in GTDB.
- Use the part AFTER the first '_' as NCBI accession (GCF_... or GCA_...).
- RS_ → use *refseq* summary; GB_ → use *genbank* summary.
- Bacteria accessions come from bac120 TSV; Archaea from ar53 TSV.
- Output is just the final HTTPS URLs (one per line), e.g. .../GCF_..._genomic.fna.gz

Stdlib only.
"""

import argparse
import csv
from pathlib import Path

FILETYPE_TO_SUFFIX = {
    "genomic_fna": "genomic.fna.gz",
    "genomic_gbff": "genomic.gbff.gz",
    "protein_faa": "protein.faa.gz",
    "gff": "genomic.gff.gz",
    "cds_from_genomic_fna": "cds_from_genomic.fna.gz",
}


def load_summary_map(path: Path) -> dict:
    """
    Return {assembly_accession -> ftp_path} from an NCBI assembly_summary.txt.
    The header line starts with '# assembly_accession\t...'; we detect it explicitly.
    """
    text = path.read_text(encoding="utf-8", errors="replace").splitlines()
    header = None
    header_idx = None

    # Find the header line (starts with '#' followed by 'assembly_accession')
    for i, line in enumerate(text):
        if line.startswith("#"):
            line2 = line[1:].lstrip()  # strip leading '#' and following spaces
            if line2.lower().startswith("assembly_accession"):
                header = line2.split("\t")
                header_idx = i
                break

    if header is None:
        raise ValueError(f"Could not find header in {path} (no '# assembly_accession' line).")

    # Build column index map
    idx = {col: j for j, col in enumerate(header)}
    if "assembly_accession" not in idx or "ftp_path" not in idx:
        raise ValueError(f"Missing required columns in {path}: need 'assembly_accession' and 'ftp_path'.")

    acc_i = idx["assembly_accession"]
    ftp_i = idx["ftp_path"]

    mapping = {}
    # Parse rows after the header; skip any additional comment lines
    for line in text[header_idx + 1:]:
        if not line or line.startswith("#"):
            continue
        row = line.split("\t")
        # guard against short/blank rows
        if len(row) <= max(acc_i, ftp_i):
            continue
        acc = row[acc_i].strip()
        ftp = row[ftp_i].strip()
        if acc and ftp:
            mapping[acc] = ftp

    return mapping


def read_gtdb_accessions(tsv_path: Path) -> list:
    """
    Read GTDB TSV and return list of (prefix, accession) where:
      prefix in {'RS','GB'}
      accession like 'GCF_...' or 'GCA_...'
    """
    out = []
    with tsv_path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        if "accession" not in reader.fieldnames:
            raise ValueError(f"{tsv_path} missing 'accession' column")
        for row in reader:
            raw = (row.get("accession") or "").strip()
            if not raw:
                continue
            parts = raw.split("_", 1)
            if len(parts) != 2:
                continue
            prefix, acc = parts[0], parts[1]
            if prefix in {"RS", "GB"} and acc:
                out.append((prefix, acc))
    # De-duplicate
    return sorted(set(out))


def make_file_url(ftp_path: str, suffix: str) -> str:
    """
    Convert ftp://... -> https://... and append '/{basename}_{suffix}'
    """
    base = ftp_path.replace("ftp://", "https://").rstrip("/")
    basename = base.rsplit("/", 1)[-1]
    return f"{base}/{basename}_{suffix}"


def main():
    ap = argparse.ArgumentParser(description="Create URL list for genome downloads.")
    ap.add_argument("--bac120", default="genome_info/bac120_metadata_filtered.tsv")
    ap.add_argument("--ar53", default="genome_info/ar53_metadata_filtered.tsv")

    ap.add_argument("--bact_refseq", default="genome_info/bacteria_refseq_assembly_summary.txt")
    ap.add_argument("--bact_genbank", default="genome_info/bacteria_genbank_assembly_summary.txt")
    ap.add_argument("--arch_refseq", default="genome_info/archaea_refseq_assembly_summary.txt")
    ap.add_argument("--arch_genbank", default="genome_info/archaea_genbank_assembly_summary.txt")

    ap.add_argument("--filetype", default="genomic_fna", choices=sorted(FILETYPE_TO_SUFFIX.keys()))
    ap.add_argument("--out_urls", default="genome_urls.txt",
                    help="Where to write the final list of URLs (one per line)")
    args = ap.parse_args()

    suffix = FILETYPE_TO_SUFFIX[args.filetype]

    # Load summaries into maps
    bact_refseq = load_summary_map(Path(args.bact_refseq))
    bact_genbank = load_summary_map(Path(args.bact_genbank))
    arch_refseq = load_summary_map(Path(args.arch_refseq))
    arch_genbank = load_summary_map(Path(args.arch_genbank))

    # Read GTDB accessions
    bact_accs = read_gtdb_accessions(Path(args.bac120))  # list of (prefix, accession)
    arch_accs = read_gtdb_accessions(Path(args.ar53))

    urls = []

    # Bacteria: RS -> RefSeq, GB -> GenBank
    for prefix, acc in bact_accs:
        ftp = bact_refseq.get(acc) if prefix == "RS" else bact_genbank.get(acc)
        if ftp:
            urls.append(make_file_url(ftp, suffix))

    # Archaea
    for prefix, acc in arch_accs:
        ftp = arch_refseq.get(acc) if prefix == "RS" else arch_genbank.get(acc)
        if ftp:
            urls.append(make_file_url(ftp, suffix))

    # De-duplicate & sort for stable output
    urls = sorted(set(urls))

    out_path = Path(args.out_urls)
    out_path.write_text("\n".join(urls) + ("\n" if urls else ""), encoding="utf-8")

    print(f"Wrote {len(urls)} URLs -> {out_path}")


if __name__ == "__main__":
    main()
