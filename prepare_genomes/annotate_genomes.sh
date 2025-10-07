#!/bin/bash

# === CONFIGURATION ===
INPUT_DIR="downloads"
OUTPUT_DIR="genomes_annotated"
THREADS=4

# === CREATE OUTPUT DIR ===
mkdir -p "$OUTPUT_DIR"

# === GET LIST OF FILES SAFELY ===
shopt -s nullglob
FILES=("$INPUT_DIR"/*.fna.gz)
shopt -u nullglob

TOTAL=${#FILES[@]}

# === CHECK FOR EMPTY INPUT ===
if [[ $TOTAL -eq 0 ]]; then
    echo "No .fna.gz files found in $INPUT_DIR. Exiting."
    exit 1
fi

# === MAIN LOOP ===
COUNT=0
for FILE in "${FILES[@]}"; do
    ((COUNT++))
    BASENAME=$(basename "$FILE" .fna.gz)
    OUTFILE="$OUTPUT_DIR/${BASENAME}.gbk"

    if [[ -f "$OUTFILE" ]]; then
        echo "Skipping $BASENAME ($COUNT/$TOTAL) — already annotated."
        continue
    fi

    echo "Annotating $BASENAME ($COUNT/$TOTAL)..."
    TMPDIR=$(mktemp -d)

    prokka "$FILE" \
        --outdir "$TMPDIR" \
        --prefix "$BASENAME" \
        --cpus "$THREADS" \
        --force > /dev/null

    mv "$TMPDIR/$BASENAME.gbk" "$OUTFILE"
    rm -r "$TMPDIR"
done

echo "All genomes processed and saved to $OUTPUT_DIR/"
