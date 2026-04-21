#!/usr/bin/env bash
set -euo pipefail

# -------- settings --------
GENOME_INFO_DIR="data/genome_info/r232"
OUT_DIR="data/genomes/r232_hq_rep_sup10"
URL_LIST="genome_urls.txt"
FILETYPE="genomic_gbff"   # genomic_fna | genomic_gbff | protein_faa | gff | cds_from_genomic_fna
PARALLEL_JOBS=8          # wget in parallel with xargs
# -------------------------

mkdir -p "$GENOME_INFO_DIR" "$OUT_DIR"

# Ensure NCBI assembly summaries exist
declare -A SUMMARY_URLS=(
  ["$GENOME_INFO_DIR/bacteria_refseq_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
  ["$GENOME_INFO_DIR/bacteria_genbank_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
  ["$GENOME_INFO_DIR/archaea_refseq_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt"
  ["$GENOME_INFO_DIR/archaea_genbank_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt"
)

for path in "${!SUMMARY_URLS[@]}"; do
  if [[ ! -s "$path" ]]; then
    echo "[DL] $path"
    wget -q -O "$path.part" "${SUMMARY_URLS[$path]}"
    mv "$path.part" "$path"
  else
    echo "[OK] $path present"
  fi
done

# Ensure GTDB metadata files exist (you provide them here)
BAC_TSV="$GENOME_INFO_DIR/bac120_metadata_filtered.tsv"
ARC_TSV="$GENOME_INFO_DIR/ar53_metadata_filtered.tsv"
if [[ ! -s "$BAC_TSV" ]] || [[ ! -s "$ARC_TSV" ]]; then
  echo "ERROR: Missing $BAC_TSV or $ARC_TSV" >&2
  exit 1
fi

# Build URL list (uses stdlib-only Python script)
python3 prepare_genomes/make_genome_url_list.py \
  --bac120 "$BAC_TSV" \
  --ar53 "$ARC_TSV" \
  --bact_refseq "$GENOME_INFO_DIR/bacteria_refseq_assembly_summary.txt" \
  --bact_genbank "$GENOME_INFO_DIR/bacteria_genbank_assembly_summary.txt" \
  --arch_refseq "$GENOME_INFO_DIR/archaea_refseq_assembly_summary.txt" \
  --arch_genbank "$GENOME_INFO_DIR/archaea_genbank_assembly_summary.txt" \
  --filetype "$FILETYPE" \
  --out_urls "$URL_LIST"

# Download all genomes into OUT_DIR.
# -c resumes partial files
# -nv quiet-ish logs
# parallel with xargs (restartable)
echo "[DL] Downloading into $OUT_DIR (parallel=$PARALLEL_JOBS)…"
xargs -a "$URL_LIST" -P "$PARALLEL_JOBS" -n 1 -I {} \
  wget -c -nv -P "$OUT_DIR" {}

echo "[DONE] All downloads attempted. Files are in: $OUT_DIR"
