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

No third-party deps (uses only csv/argparse/pathlib).
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
    Ignores lines starting with '#'.
    """
    with path.open("r", encoding="utf-8", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = None
        mapping = {}
        for row in reader:
            if not row or row[0].startswith("#"):
                # track header if present after '#'
                if row and row[0].startswith("# ") and header is None:
                    header = row[0][2:].split("\t")
                continue

            # if proper header not captured via "# ", assume first non-comment row is header
            if header is None:
                header = row
                continue

            # build index map once
            if isinstance(header, list):
                idx = {col: i for i, col in enumerate(header)}
                header = (header, idx)  # cache indices
            header_cols, idx = header

            acc = row[idx["assembly_accession"]]
            ftp = row[idx["ftp_path"]]
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
    # drop dups
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
        if prefix == "RS":
            ftp = bact_refseq.get(acc)
        else:  # GB
            ftp = bact_genbank.get(acc)
        if ftp:
            urls.append(make_file_url(ftp, suffix))

    # Archaea
    for prefix, acc in arch_accs:
        if prefix == "RS":
            ftp = arch_refseq.get(acc)
        else:
            ftp = arch_genbank.get(acc)
        if ftp:
            urls.append(make_file_url(ftp, suffix))

    # De-duplicate & sort for stable output
    urls = sorted(set(urls))

    out_path = Path(args.out_urls)
    out_path.write_text("\n".join(urls) + ("\n" if urls else ""), encoding="utf-8")

    print(f"Wrote {len(urls)} URLs -> {out_path}")

if __name__ == "__main__":
    main()
