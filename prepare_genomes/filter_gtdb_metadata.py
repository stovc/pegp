"""
Filter and Supplement GTDB Metadata by Phylum.

This script processes GTDB metadata tables (e.g. bac120, ar53) and:

1. Keeps genomes that are both:
       - mimag_high_quality == 't'
       - gtdb_representative == 't'

2. Ensures each phylum has at least TARGET_HQ_REP_MIN genomes.
   - Always retrieves additional representative genomes (if needed),
     ranked by:
         quality_score = checkm_completeness - 5 × checkm_contamination
     and filtered by QSCORE_MIN.
   - Optionally retrieves non-representative genomes
     (ENABLE_NONREP_RETRIEVAL).

Outputs (per input file):
    *_filtered.tsv            Final genome set
    *_phylum_report.xlsx      Per-phylum summary containing:
        phylum
        total
        representative
        high_quality
        high_quality_representative
        passed_filtering
        retrieved

Phylum is extracted from 'gtdb_taxonomy' (p__...).
"""

import pandas as pd
from pathlib import Path

# === CONFIGURATION ===
INPUT_FILES = [
    'data/genome_info/r232/bac120_metadata.tsv',
    'data/genome_info/r232/ar53_metadata.tsv',
]

TARGET_HQ_REP_MIN = 10       # target number of genomes per phylum
QSCORE_CONTAM_MULT = 5       # quality_score = completeness - 5 × contamination
QSCORE_MIN = 50              # minimum quality_score to be retrieved
ENABLE_NONREP_RETRIEVAL = False  # optional second retrieval step


def process_file(input_path_str: str) -> None:
    input_path = Path(input_path_str)
    df = pd.read_csv(input_path, sep='\t')

    # Ensure numeric for score computation
    for col in ['checkm_completeness', 'checkm_contamination']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # Extract phylum
    df['phylum'] = df['gtdb_taxonomy'].str.extract(r'p__([^;]+)')
    df['phylum'] = df['phylum'].fillna('Unknown')

    # === 1️⃣ Initial filtering: HQ + representative ===
    filtered_df = df[
        (df['mimag_high_quality'] == 't') &
        (df['gtdb_representative'] == 't')
    ].copy()

    # === 2️⃣ Determine phyla needing supplementation ===
    all_phyla = df['phylum'].dropna().unique()
    hq_rep_counts = (
        filtered_df.groupby('phylum').size().reset_index(name='hq_rep_count')
    )
    phyla_in_counts = set(hq_rep_counts['phylum'])
    supplement_needed = []

    for phylum in all_phyla:
        if phylum not in phyla_in_counts:
            supplement_needed.append((phylum, 0))

    supplement_needed += list(
        hq_rep_counts[hq_rep_counts['hq_rep_count'] < TARGET_HQ_REP_MIN]
        .itertuples(index=False, name=None)
    )

    # === 3️⃣ Retrieval steps ===
    additional_rows = []
    retrieved_info = []

    for phylum, hq_rep_count in supplement_needed:
        n_needed = TARGET_HQ_REP_MIN - hq_rep_count
        if n_needed <= 0:
            continue

        # Step 1: retrieve from representative genomes (not HQ)
        rep_candidates = df[
            (df['phylum'] == phylum)
            & (df['gtdb_representative'] == 't')
            & ~(df['accession'].isin(filtered_df['accession']))
        ].copy()

        rep_candidates['quality_score'] = (
            rep_candidates['checkm_completeness']
            - QSCORE_CONTAM_MULT * rep_candidates['checkm_contamination']
        )
        rep_candidates = rep_candidates[rep_candidates['quality_score'] > QSCORE_MIN]

        top_rep_candidates = rep_candidates.nlargest(n_needed, 'quality_score')
        retrieved = len(top_rep_candidates)
        n_remaining = n_needed - retrieved

        retrieved_total = retrieved
        retrieved_rows = [top_rep_candidates]

        # Step 2: optionally retrieve from non-representatives
        if ENABLE_NONREP_RETRIEVAL and n_remaining > 0:
            nonrep_candidates = df[
                (df['phylum'] == phylum)
                & (df['gtdb_representative'] != 't')
                & ~(df['accession'].isin(filtered_df['accession']))
            ].copy()

            nonrep_candidates['quality_score'] = (
                nonrep_candidates['checkm_completeness']
                - QSCORE_CONTAM_MULT * nonrep_candidates['checkm_contamination']
            )
            nonrep_candidates = nonrep_candidates[nonrep_candidates['quality_score'] > QSCORE_MIN]

            top_nonrep_candidates = nonrep_candidates.nlargest(n_remaining, 'quality_score')
            retrieved_rows.append(top_nonrep_candidates)
            retrieved_total += len(top_nonrep_candidates)

        # Combine both retrieval sources
        if retrieved_total > 0:
            additional_rows.append(pd.concat(retrieved_rows, ignore_index=True))
        retrieved_info.append((phylum, retrieved_total))

    # === 4️⃣ Combine final DataFrame ===
    if additional_rows:
        additional_df = pd.concat(additional_rows, ignore_index=True)
        result_df = pd.concat([filtered_df, additional_df], ignore_index=True)
    else:
        result_df = filtered_df.copy()

    # === 5️⃣ Save filtered TSV ===
    output_tsv = input_path.with_name(input_path.stem + '_filtered.tsv')
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"[OK] Filtered file saved: {output_tsv}")

    # === 6️⃣ Build Excel report ===
    retrieved_map = {p: n for p, n in retrieved_info}
    report_rows = []

    for ph in all_phyla:
        total = df[df['phylum'] == ph].shape[0]
        representative = df[(df['phylum'] == ph) & (df['gtdb_representative'] == 't')].shape[0]
        high_quality = df[(df['phylum'] == ph) & (df['mimag_high_quality'] == 't')].shape[0]
        high_quality_representative = df[
            (df['phylum'] == ph)
            & (df['mimag_high_quality'] == 't')
            & (df['gtdb_representative'] == 't')
        ].shape[0]
        # Passed filtering now includes all genomes in final result
        passed_filtering = result_df[result_df['phylum'] == ph].shape[0]
        retrieved = retrieved_map.get(ph, 0)

        report_rows.append({
            'phylum': ph,
            'total': total,
            'representative': representative,
            'high_quality': high_quality,
            'high_quality_representative': high_quality_representative,
            'passed_filtering': passed_filtering,
            'retrieved': retrieved,
        })

    report_df = pd.DataFrame(report_rows).sort_values('phylum').reset_index(drop=True)

    output_xlsx = input_path.with_name(input_path.stem + '_phylum_report.xlsx')
    report_df.to_excel(output_xlsx, index=False)
    print(f"[OK] Excel report saved: {output_xlsx}")


def main():
    for path in INPUT_FILES:
        print(f"\nProcessing: {path}")
        process_file(path)


if __name__ == "__main__":
    main()
