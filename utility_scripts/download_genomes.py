#!/usr/bin/env python3
"""
Download genomes listed in GTDB bac120/ar53 metadata, using NCBI assembly_summary files.

Rules:
- In the GTDB .tsv, 'accession' values look like RS_GCF_... or GB_GCA_...
- Use the part AFTER the first '_' as the true accession (e.g., GCF_..., GCA_...).
- If prefix is RS_ → use REFSEQ assembly_summary; if GB_ → use GENBANK assembly_summary.
- Bacteria accessions come from bac120 TSV; Archaea accessions come from ar53 TSV.

Outputs:
- All genomes saved into: out_dir/ (single folder)
- Report CSV: out_dir/download_report.csv
- All metadata (GTDB TSV + assembly summaries) are under: genome_info/
"""

import argparse
import concurrent.futures as cf
import csv
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import requests


# =========================
# Constants
# =========================

FILETYPE_TO_SUFFIX = {
    "genomic_fna": "genomic.fna.gz",
    "genomic_gbff": "genomic.gbff.gz",
    "protein_faa": "protein.faa.gz",
    "gff": "genomic.gff.gz",
    "cds_from_genomic_fna": "cds_from_genomic.fna.gz",
}

SUMMARY_URLS = {
    "bacteria_refseq": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt",
    "bacteria_genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt",
    "archaea_refseq": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt",
    "archaea_genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt",
}

# Tunables (can be overridden by CLI)
max_workers = 8
retries = 3
timeout = 60  # seconds
backoff = 2.0  # exponential backoff base


# =========================
# Helpers
# =========================

def ensure_dir(path: Path):
    """Create directory (and parents) if it doesn't exist."""
    path.mkdir(parents=True, exist_ok=True)


def ensure_summary_file(path: Path, key: str):
    """Download NCBI assembly_summary.txt if it does not already exist (non-empty)."""
    if path.exists() and path.stat().st_size > 0:
        print(f"[OK] Found summary: {path}")
        return

    url = SUMMARY_URLS[key]
    print(f"[DL] Downloading {key} from {url}")
    ensure_dir(path.parent)
    with requests.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    print(f"[DONE] Saved {path}")


def load_summary(path: Path) -> pd.DataFrame:
    """Load NCBI assembly_summary .txt (tab-delimited, comment lines start with '#')."""
    df = pd.read_csv(
        path,
        sep='\t',
        comment='#',
        dtype=str,
        low_memory=False
    )
    needed = ['assembly_accession', 'ftp_path', 'version_status', 'asm_name']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")
    return df[needed].dropna(subset=['assembly_accession', 'ftp_path'])


def read_accessions(tsv_path: Path) -> pd.DataFrame:
    """
    Read a GTDB metadata TSV and return columns:
    ['accession_raw','prefix','accession'].

    'accession' in GTDB looks like 'RS_GCF_...' or 'GB_GCA_...'.
    We split on the first underscore: prefix = 'RS'|'GB', accession = 'GCF_...|GCA_...'.
    """
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)
    if 'accession' not in df.columns:
        raise ValueError(f"File {tsv_path} has no 'accession' column.")
    acc = df['accession'].astype(str)
    parts = acc.str.split('_', n=1, expand=True)
    prefix = parts[0]
    acc2 = parts[1]
    out = pd.DataFrame({
        'accession_raw': acc,
        'prefix': prefix,
        'accession': acc2
    })
    out = out[out['prefix'].isin(['RS', 'GB']) & out['accession'].notna()].drop_duplicates()
    return out


def make_url(ftp_path: str, filetype: str) -> str:
    """
    Build HTTPS URL to a specific file under an assembly ftp_path.
    ftp_path like: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/.../GCF_XXXX.Y_ASM_name
    We convert ftp->https and append '{basename}_{suffix}'.
    """
    suffix = FILETYPE_TO_SUFFIX[filetype]
    base = ftp_path.replace("ftp://", "https://").rstrip('/')
    basename = base.split('/')[-1]
    return f"{base}/{basename}_{suffix}"


def download_file(url: str, dest: Path) -> Tuple[bool, str]:
    """Download with retries. Returns (success, message)."""
    for attempt in range(1, retries + 1):
        try:
            with requests.get(url, stream=True, timeout=timeout) as r:
                if r.status_code != 200:
                    raise RuntimeError(f"HTTP {r.status_code}")
                ensure_dir(dest.parent)
                with open(dest, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
            return True, "ok"
        except Exception as e:
            if attempt == retries:
                return False, f"failed: {e}"
            time.sleep(backoff ** attempt)
    return False, "unknown error"


# =========================
# Core logic
# =========================

def build_lookup(
    bacteria_refseq_summary: Path,
    bacteria_genbank_summary: Path,
    archaea_refseq_summary: Path,
    archaea_genbank_summary: Path,
) -> Dict[Tuple[str, str], pd.DataFrame]:
    """Prepare summary lookups keyed by (kingdom, source)."""
    return {
        ('bacteria', 'refseq'): load_summary(bacteria_refseq_summary),
        ('bacteria', 'genbank'): load_summary(bacteria_genbank_summary),
        ('archaea', 'refseq'): load_summary(archaea_refseq_summary),
        ('archaea', 'genbank'): load_summary(archaea_genbank_summary),
    }


def resolve_paths(
    accessions_df: pd.DataFrame,
    kingdom: str,
    lookup: Dict[Tuple[str, str], pd.DataFrame],
) -> pd.DataFrame:
    rs = accessions_df[accessions_df['prefix'] == 'RS'].copy()
    gb = accessions_df[accessions_df['prefix'] == 'GB'].copy()

    out_parts = []

    if not rs.empty:
        df_sum = lookup[(kingdom, 'refseq')]
        merged = rs.merge(df_sum, left_on='accession', right_on='assembly_accession', how='left')
        merged['source'] = 'refseq'
        out_parts.append(merged)

    if not gb.empty:
        df_sum = lookup[(kingdom, 'genbank')]
        merged = gb.merge(df_sum, left_on='accession', right_on='assembly_accession', how='left')
        merged['source'] = 'genbank'
        out_parts.append(merged)

    if out_parts:
        out = pd.concat(out_parts, ignore_index=True)
    else:
        out = pd.DataFrame(columns=[
            'accession_raw', 'prefix', 'accession',
            'assembly_accession', 'ftp_path', 'version_status',
            'asm_name', 'source'
        ])

    out['kingdom'] = kingdom
    return out


def plan_downloads(
    resolved_df: pd.DataFrame,
    filetype: str,
    out_dir: Path,
) -> List[Dict]:
    """Produce a list of planned downloads with url and destination (single folder)."""
    plans = []
    for _, row in resolved_df.iterrows():
        record = {
            'kingdom': row.get('kingdom', ''),
            'source': row.get('source', ''),
            'accession_raw': row.get('accession_raw', ''),
            'accession': row.get('accession', ''),
            'ftp_path': row.get('ftp_path', None),
            'url': None,
            'dest': None,
        }

        ftp = record['ftp_path']
        if pd.isna(ftp) or not isinstance(ftp, str) or not ftp.strip():
            plans.append(record)
            continue

        try:
            url = make_url(ftp, filetype)
            # Save everything into a single folder
            dest = out_dir / f"{url.split('/')[-1]}"
            record['url'] = url
            record['dest'] = str(dest)
        except Exception:
            pass

        plans.append(record)
    return plans


def run_downloads(plans: List[Dict]) -> List[Dict]:
    """
    Execute downloads in parallel; skip files that already exist and are non-empty.
    Returns list of results with status fields: ok | skipped | missing_ftp | error
    """
    results = []

    def _task(plan: Dict):
        url = plan['url']
        dest = plan['dest']
        if not url or not dest:
            return {**plan, 'status': 'missing_ftp', 'message': 'No ftp_path/url'}

        dest_path = Path(dest)
        if dest_path.exists() and dest_path.stat().st_size > 0:
            return {**plan, 'status': 'skipped', 'message': 'already exists'}

        ok, msg = download_file(url, dest_path)
        return {**plan, 'status': 'ok' if ok else 'error', 'message': msg}

    with cf.ThreadPoolExecutor(max_workers=max_workers) as ex:
        for res in ex.map(_task, plans):
            results.append(res)
    return results


def write_report(rows: List[Dict], path: Path):
    """Save a CSV report of actions taken."""
    ensure_dir(path.parent)
    cols = ['kingdom', 'source', 'accession_raw', 'accession', 'ftp_path', 'url', 'dest', 'status', 'message']
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=cols)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k, '') for k in cols})


# =========================
# CLI
# =========================

def main():
    parser = argparse.ArgumentParser(description="Download NCBI genomes listed in GTDB bac120/ar53 metadata.")

    # All metadata defaults under genome_info/
    parser.add_argument("--bac120", default="genome_info/bac120_metadata_filtered.tsv",
                        help="Path to GTDB bacteria TSV (bac120)")
    parser.add_argument("--ar53", default="genome_info/ar53_metadata_filtered.tsv",
                        help="Path to GTDB archaea TSV (ar53)")

    parser.add_argument("--bact_refseq", default="genome_info/bacteria_refseq_assembly_summary.txt",
                        help="Bacteria RefSeq assembly_summary")
    parser.add_argument("--bact_genbank", default="genome_info/bacteria_genbank_assembly_summary.txt",
                        help="Bacteria GenBank assembly_summary")
    parser.add_argument("--arch_refseq", default="genome_info/archaea_refseq_assembly_summary.txt",
                        help="Archaea RefSeq assembly_summary")
    parser.add_argument("--arch_genbank", default="genome_info/archaea_genbank_assembly_summary.txt",
                        help="Archaea GenBank assembly_summary")

    parser.add_argument("--out_dir", default="downloads",
                        help="Output directory for downloaded genomes (single folder)")
    parser.add_argument("--filetype", default="genomic_fna",
                        choices=sorted(FILETYPE_TO_SUFFIX.keys()),
                        help="Which file to download (default: genomic_fna → *_genomic.fna.gz)")
    parser.add_argument("--workers", type=int, default=max_workers, help="Parallel downloads")
    parser.add_argument("--retries", type=int, default=retries, help="Download retries")

    args = parser.parse_args()

    # Apply CLI overrides
    global max_workers, retries
    max_workers = args.workers
    retries = args.retries

    out_dir = Path(args.out_dir)
    filetype = args.filetype
    ensure_dir(out_dir)

    # Ensure genome_info/ exists (based on defaults) to place summaries if needed
    genome_info_dir = Path(args.bact_refseq).parent
    ensure_dir(genome_info_dir)

    # Ensure assembly summaries exist (download if missing)
    ensure_summary_file(Path(args.bact_refseq), "bacteria_refseq")
    ensure_summary_file(Path(args.bact_genbank), "bacteria_genbank")
    ensure_summary_file(Path(args.arch_refseq), "archaea_refseq")
    ensure_summary_file(Path(args.arch_genbank), "archaea_genbank")

    # Build lookup tables
    lookup = build_lookup(
        Path(args.bact_refseq),
        Path(args.bact_genbank),
        Path(args.arch_refseq),
        Path(args.arch_genbank)
    )

    # Read GTDB accessions
    bac_df = read_accessions(Path(args.bac120))
    arc_df = read_accessions(Path(args.ar53))

    # Resolve download paths
    bac_resolved = resolve_paths(bac_df, 'bacteria', lookup)
    arc_resolved = resolve_paths(arc_df, 'archaea', lookup)

    # Plan & download (all into a single folder)
    plans = plan_downloads(bac_resolved, filetype, out_dir) + plan_downloads(arc_resolved, filetype, out_dir)
    if not plans:
        print("No accessions to process.", file=sys.stderr)
        sys.exit(1)

    results = run_downloads(plans)

    # Write report
    report_path = out_dir / "download_report.csv"
    write_report(results, report_path)

    # Summary to stdout
    ok = sum(1 for r in results if r.get('status') == 'ok')
    skipped = sum(1 for r in results if r.get('status') == 'skipped')
    missing = sum(1 for r in results if r.get('status') == 'missing_ftp')
    err = sum(1 for r in results if r.get('status') == 'error')
    print(f"Done. OK={ok}, skipped={skipped}, missing_ftp={missing}, error={err}")
    print(f"Report: {report_path.resolve()}")


if __name__ == "__main__":
    main()
