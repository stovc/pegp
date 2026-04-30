"""
Filter GTDB Metadata by Configurable Genome-Quality Criteria.

This script processes GTDB metadata tables and keeps genomes according to
configurable filters:
    - mimag_high_quality == 't'
    - ncbi_assembly_level == 'Complete Genome'
    - gtdb_representative == 't'

Outputs per input file:
    *_filtered.tsv
    *_phylum_report.xlsx

The Excel report summarizes genome counts per phylum and in total.
"""

from pathlib import Path

import pandas as pd


# === CONFIGURATION ===
INPUT_FILES = [
    'data/genome_info/r232/bac120_metadata.tsv',
    'data/genome_info/r232/ar53_metadata.tsv',
]

FILTER_MIMAG_HIGH_QUALITY = True
FILTER_COMPLETE_GENOME = True
FILTER_GTDB_REPRESENTATIVE = False


def process_file(input_path_str: str) -> None:
    input_path = Path(input_path_str)
    df = pd.read_csv(input_path, sep='\t')

    # Extract phylum
    df['phylum'] = df['gtdb_taxonomy'].str.extract(r'p__([^;]+)')
    df['phylum'] = df['phylum'].fillna('Unknown')

    # === Apply configurable filters ===
    filter_mask = pd.Series(True, index=df.index)

    if FILTER_MIMAG_HIGH_QUALITY:
        filter_mask &= df['mimag_high_quality'] == 't'

    if FILTER_COMPLETE_GENOME:
        filter_mask &= df['ncbi_assembly_level'] == 'Complete Genome'

    if FILTER_GTDB_REPRESENTATIVE:
        filter_mask &= df['gtdb_representative'] == 't'

    result_df = df[filter_mask].copy()

    # === Save filtered TSV ===
    output_tsv = input_path.with_name(input_path.stem + '_filtered.tsv')
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"[OK] Filtered file saved: {output_tsv}")

    # === Build report ===
    report_rows = []

    for phylum, phylum_df in df.groupby('phylum', dropna=False):
        result_phylum_df = result_df[result_df['phylum'] == phylum]

        total = phylum_df.shape[0]
        mimag_high_quality = (phylum_df['mimag_high_quality'] == 't').sum()
        complete = (phylum_df['ncbi_assembly_level'] == 'Complete Genome').sum()
        complete_high_quality = (
            (phylum_df['mimag_high_quality'] == 't')
            & (phylum_df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum()
        passed_filtering = result_phylum_df.shape[0]

        report_rows.append({
            'phylum': phylum,
            'total': total,
            'mimag_high_quality': mimag_high_quality,
            'complete': complete,
            'complete_high_quality': complete_high_quality,
            'passed_filtering': passed_filtering,
        })

    # Add total row
    report_rows.append({
        'phylum': 'TOTAL',
        'total': df.shape[0],
        'mimag_high_quality': (df['mimag_high_quality'] == 't').sum(),
        'complete': (df['ncbi_assembly_level'] == 'Complete Genome').sum(),
        'complete_high_quality': (
            (df['mimag_high_quality'] == 't')
            & (df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum(),
        'passed_filtering': result_df.shape[0],
    })

    report_df = (
        pd.DataFrame(report_rows)
        .sort_values('total', ascending=False)
        .reset_index(drop=True)
    )

    output_xlsx = input_path.with_name(input_path.stem + '_phylum_report.xlsx')
    report_df.to_excel(output_xlsx, index=False)
    print(f"[OK] Excel report saved: {output_xlsx}")


def main() -> None:
    for path in INPUT_FILES:
        print(f"\nProcessing: {path}")
        process_file(path)


if __name__ == "__main__":
    main()