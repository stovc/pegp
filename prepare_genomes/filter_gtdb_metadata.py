"""
Filter GTDB Metadata by Configurable Genome-Quality Criteria.

This script processes GTDB metadata tables and keeps genomes according to
configurable filters:
    - mimag_high_quality == 't'
    - ncbi_assembly_level == 'Complete Genome'
    - gtdb_representative == 't'
    - gtdb_type_species_of_genus == 't'
    - database inferred from accession prefix: RS or GB

Optionally, after all filters are applied, the script keeps only N genomes
per selected taxonomic rank, prioritizing type species, representatives,
RefSeq records, complete genomes, and higher quality_score.

Taxa represented by fewer than a configurable number of eligible genomes can
also be removed before the per-taxon limit is applied.

Outputs per input file:
    *_filtered.tsv
    *_phylum_report.xlsx
"""

from pathlib import Path

import pandas as pd


# === CONFIGURATION ===
INPUT_FILES = [
    'data/genome_info/r232/bac120_metadata.tsv',
    'data/genome_info/r232/ar53_metadata.tsv',
]

FILTER_MIMAG_HIGH_QUALITY = False
FILTER_COMPLETE_GENOME = True
FILTER_GTDB_REPRESENTATIVE = True
FILTER_GTDB_TYPE_SPECIES_OF_GENUS = True

# Options:
#   None        -> no database filtering
#   'RS'        -> keep only RefSeq
#   'GB'        -> keep only GenBank
#   ['RS','GB'] -> keep both; equivalent to no database filtering if only these exist
DATABASE_FILTER = 'RS'

# Applied after the filters above, before limiting genomes per taxon
ENABLE_MIN_GENOMES_PER_TAXON = True
MIN_GENOMES_TAXONOMIC_RANK = 'phylum'
MIN_GENOMES_PER_TAXON = 3

# Applied last, after the filters above
ENABLE_LIMIT_PER_TAXON = True
TAXONOMIC_RANK = 'family'
GENOMES_PER_TAXON = 1

QSCORE_CONTAM_MULT = 5


RANK_PREFIXES = {
    'domain': 'd__',
    'phylum': 'p__',
    'class': 'c__',
    'order': 'o__',
    'family': 'f__',
    'genus': 'g__',
    'species': 's__',
}


def extract_taxon(taxonomy: pd.Series, rank: str) -> pd.Series:
    if rank not in RANK_PREFIXES:
        valid_ranks = ', '.join(RANK_PREFIXES)
        raise ValueError(f"Invalid rank '{rank}'. Valid ranks: {valid_ranks}")

    prefix = RANK_PREFIXES[rank]
    return taxonomy.str.extract(fr'{prefix}([^;]+)')[0].fillna('Unknown')


def normalize_database_filter(database_filter):
    if database_filter is None:
        return None

    if isinstance(database_filter, str):
        return [database_filter]

    return list(database_filter)


def add_helper_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    df['phylum'] = extract_taxon(df['gtdb_taxonomy'], 'phylum')
    df['selected_taxon'] = extract_taxon(df['gtdb_taxonomy'], TAXONOMIC_RANK)

    # Database is inferred from accession prefix: RS_..., GB_...
    df['database'] = df['accession'].astype(str).str[:2]

    for col in ['checkm_completeness', 'checkm_contamination']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df['quality_score'] = (
        df['checkm_completeness']
        - QSCORE_CONTAM_MULT * df['checkm_contamination']
    )

    return df


def apply_basic_filters(df: pd.DataFrame) -> pd.DataFrame:
    filter_mask = pd.Series(True, index=df.index)

    if FILTER_MIMAG_HIGH_QUALITY:
        filter_mask &= df['mimag_high_quality'] == 't'

    if FILTER_COMPLETE_GENOME:
        filter_mask &= df['ncbi_assembly_level'] == 'Complete Genome'

    if FILTER_GTDB_REPRESENTATIVE:
        filter_mask &= df['gtdb_representative'] == 't'

    if FILTER_GTDB_TYPE_SPECIES_OF_GENUS:
        filter_mask &= df['gtdb_type_species_of_genus'] == 't'

    database_filter = normalize_database_filter(DATABASE_FILTER)
    if database_filter is not None:
        filter_mask &= df['database'].isin(database_filter)

    return df[filter_mask].copy()


def filter_taxa_by_min_genomes(df: pd.DataFrame) -> pd.DataFrame:
    if not ENABLE_MIN_GENOMES_PER_TAXON:
        return df.copy()

    if MIN_GENOMES_PER_TAXON < 1:
        raise ValueError('MIN_GENOMES_PER_TAXON must be at least 1')

    taxa = extract_taxon(
        df['gtdb_taxonomy'],
        MIN_GENOMES_TAXONOMIC_RANK,
    )
    taxon_counts = taxa.groupby(taxa, dropna=False).transform('size')

    return df[taxon_counts >= MIN_GENOMES_PER_TAXON].copy()


def limit_genomes_per_taxon(df: pd.DataFrame) -> pd.DataFrame:
    if not ENABLE_LIMIT_PER_TAXON:
        return df.copy()

    df = df.copy()

    # Priority columns.
    # Note: user request repeated gtdb_type_species_of_genus twice.
    # I interpret the second priority as gtdb_representative.
    df['_priority_type_species'] = df['gtdb_type_species_of_genus'] == 't'
    df['_priority_representative'] = df['gtdb_representative'] == 't'
    df['_priority_rs'] = df['database'] == 'RS'
    df['_priority_complete'] = df['ncbi_assembly_level'] == 'Complete Genome'

    df = df.sort_values(
        by=[
            'selected_taxon',
            '_priority_type_species',
            '_priority_representative',
            '_priority_rs',
            '_priority_complete',
            'quality_score',
        ],
        ascending=[
            True,
            False,
            False,
            False,
            False,
            False,
        ],
    )

    limited_df = (
        df.groupby('selected_taxon', group_keys=False, dropna=False)
        .head(GENOMES_PER_TAXON)
        .copy()
    )

    limited_df = limited_df.drop(
        columns=[
            '_priority_type_species',
            '_priority_representative',
            '_priority_rs',
            '_priority_complete',
        ]
    )

    return limited_df


def build_report(df: pd.DataFrame, result_df: pd.DataFrame) -> pd.DataFrame:
    report_rows = []

    for phylum, phylum_df in df.groupby('phylum', dropna=False):
        result_phylum_df = result_df[result_df['phylum'] == phylum]

        total = phylum_df.shape[0]

        mimag_high_quality = (phylum_df['mimag_high_quality'] == 't').sum()

        complete = (
            phylum_df['ncbi_assembly_level'] == 'Complete Genome'
        ).sum()

        complete_high_quality = (
            (phylum_df['mimag_high_quality'] == 't')
            & (phylum_df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum()

        mag = (phylum_df['ncbi_genome_category'] != 'none').sum()

        gb = (phylum_df['database'] == 'GB').sum()
        rs = (phylum_df['database'] == 'RS').sum()

        complete_rs = (
            (phylum_df['database'] == 'RS')
            & (phylum_df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum()

        complete_high_quality_rs = (
            (phylum_df['database'] == 'RS')
            & (phylum_df['mimag_high_quality'] == 't')
            & (phylum_df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum()

        passed_filtering = result_phylum_df.shape[0]

        passed_filtering_mag = (
            result_phylum_df['ncbi_genome_category'] != 'none'
        ).sum()

        passed_filtering_gb = (
            result_phylum_df['database'] == 'GB'
        ).sum()

        passed_filtering_not_hq = (
            result_phylum_df['mimag_high_quality'] != 't'
        ).sum()

        passed_filtering_not_complete = (
            result_phylum_df['ncbi_assembly_level'] != 'Complete Genome'
        ).sum()

        report_rows.append({
            'phylum': phylum,
            'total': total,
            'mimag_high_quality': mimag_high_quality,
            'complete': complete,
            'complete_high_quality': complete_high_quality,
            'mag': mag,
            'GB': gb,
            'RS': rs,
            'complete_RS': complete_rs,
            'complete_high_quality_RS': complete_high_quality_rs,
            'passed_filtering': passed_filtering,
            'passed_filtering_mag': passed_filtering_mag,
            'passed_filtering_GB': passed_filtering_gb,
            'passed_filtering_not_hq': passed_filtering_not_hq,
            'passed_filtering_not_complete': passed_filtering_not_complete,
        })

    total_row = {
        'phylum': 'TOTAL',
        'total': df.shape[0],
        'mimag_high_quality': (df['mimag_high_quality'] == 't').sum(),
        'complete': (df['ncbi_assembly_level'] == 'Complete Genome').sum(),
        'complete_high_quality': (
            (df['mimag_high_quality'] == 't')
            & (df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum(),
        'mag': (df['ncbi_genome_category'] != 'none').sum(),
        'GB': (df['database'] == 'GB').sum(),
        'RS': (df['database'] == 'RS').sum(),
        'complete_RS': (
            (df['database'] == 'RS')
            & (df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum(),
        'complete_high_quality_RS': (
            (df['database'] == 'RS')
            & (df['mimag_high_quality'] == 't')
            & (df['ncbi_assembly_level'] == 'Complete Genome')
        ).sum(),
        'passed_filtering': result_df.shape[0],
        'passed_filtering_mag': (
            result_df['ncbi_genome_category'] != 'none'
        ).sum(),
        'passed_filtering_GB': (result_df['database'] == 'GB').sum(),
        'passed_filtering_not_hq': (
            result_df['mimag_high_quality'] != 't'
        ).sum(),
        'passed_filtering_not_complete': (
            result_df['ncbi_assembly_level'] != 'Complete Genome'
        ).sum(),
    }

    report_df = pd.DataFrame(report_rows)
    report_df = report_df.sort_values('total', ascending=False).reset_index(drop=True)

    report_df = pd.concat(
        [pd.DataFrame([total_row]), report_df],
        ignore_index=True,
    )

    return report_df


def process_file(input_path_str: str) -> None:
    input_path = Path(input_path_str)

    df = pd.read_csv(input_path, sep='\t')
    df = add_helper_columns(df)

    filtered_df = apply_basic_filters(df)
    filtered_df = filter_taxa_by_min_genomes(filtered_df)
    result_df = limit_genomes_per_taxon(filtered_df)

    output_tsv = input_path.with_name(input_path.stem + '_filtered.tsv')
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"[OK] Filtered file saved: {output_tsv}")

    report_df = build_report(df, result_df)

    output_xlsx = input_path.with_name(input_path.stem + '_phylum_report.xlsx')
    report_df.to_excel(output_xlsx, index=False)
    print(f"[OK] Excel report saved: {output_xlsx}")


def main() -> None:
    for path in INPUT_FILES:
        print(f"\nProcessing: {path}")
        process_file(path)


if __name__ == "__main__":
    main()
