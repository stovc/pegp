#!/usr/bin/env bash
# Combined download script: NCBI (assembly_summaries) + GTDB metadata
set -euo pipefail

OUTDIR="genomes_info"
LOGFILE="${OUTDIR}/download.log"
mkdir -p "$OUTDIR"
exec > >(tee -a "$LOGFILE") 2>&1

echo "=== Start downloads ==="
echo "Timestamp: $(date)"
echo "Output dir: $OUTDIR"
echo "---------------------------"

declare -A FILES=(
  # NCBI assembly summaries
  ["bacteria_refseq_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
  ["bacteria_genbank_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt"
  ["archaea_refseq_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt"
  ["archaea_genbank_assembly_summary.txt"]="https://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt"
  # GTDB metadata (gzipped)
  ["bac120_metadata.tsv.gz"]="https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tsv.gz"
  ["ar53_metadata.tsv.gz"]="https://data.gtdb.ecogenomic.org/releases/latest/ar53_metadata.tsv.gz"
)

for FNAME in "${!FILES[@]}"; do
    URL="${FILES[$FNAME]}"
    DEST="$OUTDIR/$FNAME"

    if [ -s "$DEST" ]; then
        echo "[skip] $FNAME exists ($(du -h "$DEST" | cut -f1))"
        continue
    fi

    echo "[download] $FNAME"
    wget -q -c --tries=3 --timeout=60 -O "$DEST" "$URL" \
      && echo "[ok] $FNAME saved ($(du -h "$DEST" | cut -f1))" \
      || { echo "[error] Failed to fetch $URL"; rm -f "$DEST"; exit 1; }
done

echo "=== Downloads done ==="
echo "You have the following files in $OUTDIR:"
ls -l "$OUTDIR"

# Optional: decompress the metadata files (if you want uncompressed .tsv)
# echo "Decompressing GTDB metadata..."
# for f in bac120_metadata.tsv.gz ar53_metadata.tsv.gz; do
#   if [ -f "$OUTDIR/$f" ]; then
#     gunzip -kf "$OUTDIR/$f"
#   fi
# done

echo "=== Finished ==="
