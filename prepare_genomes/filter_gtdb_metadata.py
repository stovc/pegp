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

Taxa represented by fewer than a configurable number of genomes can also be
removed from the final limited result.

Outputs per input file:
    *_filtered.tsv
    *_phylum_report.xlsx
"""

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_INPUT_FILES = [
    'data/genome_info/r232/bac120_metadata.tsv',
    'data/genome_info/r232/ar53_metadata.tsv',
]

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


def add_helper_columns(
    df: pd.DataFrame,
    taxonomic_rank: str,
    qscore_contam_mult: float,
) -> pd.DataFrame:
    df = df.copy()

    df['phylum'] = extract_taxon(df['gtdb_taxonomy'], 'phylum')
    df['selected_taxon'] = extract_taxon(df['gtdb_taxonomy'], taxonomic_rank)

    # Database is inferred from accession prefix: RS_..., GB_...
    df['database'] = df['accession'].astype(str).str[:2]

    for col in ['checkm_completeness', 'checkm_contamination']:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df['quality_score'] = (
        df['checkm_completeness']
        - qscore_contam_mult * df['checkm_contamination']
    )

    return df


def apply_basic_filters(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    filter_mask = pd.Series(True, index=df.index)

    if args.filter_mimag_high_quality:
        filter_mask &= df['mimag_high_quality'] == 't'

    if args.filter_complete_genome:
        filter_mask &= df['ncbi_assembly_level'] == 'Complete Genome'

    if args.filter_gtdb_representative:
        filter_mask &= df['gtdb_representative'] == 't'

    if args.filter_gtdb_type_species_of_genus:
        filter_mask &= df['gtdb_type_species_of_genus'] == 't'

    if args.database_filter:
        filter_mask &= df['database'].isin(args.database_filter)

    return df[filter_mask].copy()


def filter_taxa_by_min_genomes(
    df: pd.DataFrame,
    enabled: bool,
    rank: str,
    minimum: int,
) -> pd.DataFrame:
    if not enabled:
        return df.copy()

    taxa = extract_taxon(df['gtdb_taxonomy'], rank)
    taxon_counts = taxa.groupby(taxa, dropna=False).transform('size')

    return df[taxon_counts >= minimum].copy()


def limit_genomes_per_taxon(
    df: pd.DataFrame,
    enabled: bool,
    genomes_per_taxon: int,
) -> pd.DataFrame:
    if not enabled:
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
        .head(genomes_per_taxon)
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


def process_file(input_path_str: str, args: argparse.Namespace) -> None:
    input_path = Path(input_path_str)

    df = pd.read_csv(input_path, sep='\t')
    df = add_helper_columns(df, args.taxonomic_rank, args.qscore_contam_mult)

    filtered_df = apply_basic_filters(df, args)
    result_df = limit_genomes_per_taxon(
        filtered_df, args.enable_limit_per_taxon, args.genomes_per_taxon
    )
    result_df = filter_taxa_by_min_genomes(
        result_df,
        args.enable_min_genomes_per_taxon,
        args.min_genomes_taxonomic_rank,
        args.min_genomes_per_taxon,
    )

    output_tsv = input_path.with_name(input_path.stem + '_filtered.tsv')
    result_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"[OK] Filtered file saved: {output_tsv}")

    report_df = build_report(df, result_df)

    output_xlsx = input_path.with_name(input_path.stem + '_phylum_report.xlsx')
    report_df.to_excel(output_xlsx, index=False)
    print(f"[OK] Excel report saved: {output_xlsx}")


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError('must be at least 1')
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'input_files', nargs='*', default=DEFAULT_INPUT_FILES,
        help='GTDB metadata TSV files (default: r232 bacterial and archaeal files)',
    )

    boolean_options = (
        ('filter-mimag-high-quality', False, 'filter for MIMAG high-quality genomes'),
        ('filter-complete-genome', True, 'filter for complete genomes'),
        ('filter-gtdb-representative', True, 'filter for GTDB representatives'),
        ('filter-gtdb-type-species-of-genus', True, 'filter for GTDB type species'),
        ('enable-min-genomes-per-taxon', True, 'remove taxa below the minimum size'),
        ('enable-limit-per-taxon', True, 'limit the number of genomes per taxon'),
    )
    for name, default, help_text in boolean_options:
        parser.add_argument(
            f'--{name}', action=argparse.BooleanOptionalAction,
            default=default, help=help_text,
        )

    parser.add_argument(
        '--database-filter', nargs='+', choices=('RS', 'GB'), default=['RS'],
        metavar='{RS,GB}',
        help='databases to retain; use --no-database-filter to disable (default: RS)',
    )
    parser.add_argument(
        '--no-database-filter', dest='database_filter', action='store_const',
        const=None, help=argparse.SUPPRESS,
    )
    parser.add_argument(
        '--min-genomes-taxonomic-rank', choices=RANK_PREFIXES, default='phylum',
    )
    parser.add_argument('--min-genomes-per-taxon', type=positive_int, default=3)
    parser.add_argument('--taxonomic-rank', choices=RANK_PREFIXES, default='family')
    parser.add_argument('--genomes-per-taxon', type=positive_int, default=1)
    parser.add_argument('--qscore-contam-mult', type=float, default=5)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    for path in args.input_files:
        print(f"\nProcessing: {path}")
        process_file(path, args)


if __name__ == "__main__":
    main()
