#!/bin/bash

# === CONFIGURATION ===
INPUT_DIR="genomes"
OUTPUT_DIR="genomes_annotated"
THREADS=4

# === CREATE OUTPUT DIR ===
mkdir -p "$OUTPUT_DIR"

# === COLLECT INPUT FILES ===
FILES=("$INPUT_DIR"/*.fna.gz)
TOTAL=${#FILES[@]}
COUNT=0

# === MAIN LOOP ===
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
