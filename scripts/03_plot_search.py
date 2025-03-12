"""Plot HMMER hits and their properties.

- input: "hits_df.csv"
- output: "search_report.pdf"

Output contains the following:
    - Information about database, query, homology search tool, and numerical characteristics of search results
    - Histogram of unique genomes
    - Histogram of unique proteins
    - Histogram of hits pairwise overlaps lengthes
    - 3d plot of identity, query_coverage, and length colored by taxon
    - Pair scatter plot of identity, query_coverage, and length with histograms at the diagonal
    - Histogram of total hits per phylum (TODO: would be nice to plot values normalized to numbers of genomes/phylum)
    - Histogram of average number of hits per genome in a phylum (TODO: would be nice to plot a tree with pattern)

    TODO: make colormap generation less hardcoded (infinite value-proof)
    TODO: fix plot_overlaps()
"""

from Bio import SeqIO, SearchIO
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.backends.backend_pdf
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import seaborn as sns

import argparse
from collections import defaultdict
import io
from pathlib import Path
import sys
import traceback

import constants as c
from plot_search_utils import *


def filter_top_column_values(df, filter_column ='phylum', assign_column='taxon', keep_first=12):
    """
    Keeps the top `keep_first` most abundant taxa from and groups the rest as 'Other'.
    Parameters:
        df : DataFrame containing a column with taxonomic ranks.
        filter_column : Column name of the taxonomic rank (default: 'phylum').
        assign_column : Column name of the assigned filtered column.
        keep_first : Number of top taxa to retain (default: 12).
    Returns: DataFrame with an added 'Taxon' column, grouping rare taxa as 'Other'.
    """
    value_counts = df[filter_column].value_counts()
    top_n_values = value_counts.nlargest(keep_first).index
    df[assign_column] = df[filter_column].apply(lambda x: x if x in top_n_values else 'Other')

    return df


def save_plots_to_files(plots):
    for i, fig in enumerate(plots):
        fig.savefig(f'plot_{i+1}.png')
        plt.close(fig)  # Free memory


def save_to_pdf(text_lines, plots, pdf_path):
    c = canvas.Canvas(pdf_path, pagesize=letter)

    # Add text report to first page
    c.setFont("Helvetica", 12)
    y_pos = 750
    for line in text_lines:
        c.drawString(100, y_pos, line)
        y_pos -= 20
    c.showPage()  # Move to next page

    # Add plots to subsequent pages
    for fig in plots:
        img_buffer = io.BytesIO()
        fig.savefig(img_buffer, format='png')
        img_buffer.seek(0)
        plt.close(fig)  # Free memory
        c.drawImage(ImageReader(img_buffer), 100, 400, width=400, height=300)
        c.showPage()

    c.save()


def generate_text_report(df_handle, hmmer_results):
    protein_faa_path = Path('databases') / database / 'protein.faa'
    database_size = count_fasta_records(protein_faa_path)
    hit_genomes_number = len(df_handle['assembly'].unique())
    hit_proteins_number = len(df_handle['lcs'].unique())
    number_of_hsps = len(df_handle.index)

    text_report = [
        f"Homology search using: {hmmer_results.program} {hmmer_results.version}",
        f"Query: {hmmer_results._id}",
        f"Database: {database}",
        "Preliminary e-value cutoff: < 0.1",
        f"Database size: {database_size} proteins",
        f"Database size: {database_size} proteins",
        f"Hit genomes: {hit_genomes_number}",
        f"Hit proteins: {hit_proteins_number}",
        f"HSPs: {number_of_hsps}"
    ]
    return text_report


def pairplot(data, columns, hue, palette):
    """Plot seaborn pairplot and return a figure.

    Parameters:
    - data: DataFrame
    - columns: List of columns from data to be plotted
    - hue: Column that defines colors of plotted dots
    - palette: Color palette for dots

    Returns:
    - fig: Matplotlib figure object
    """
    g = sns.PairGrid(data[columns], diag_sharey=False, hue=hue,
                     height=3, aspect=1.5, palette=palette, corner=True)

    g.map_lower(sns.scatterplot, edgecolor="k", linewidth=0.2, s=3)
    g.map_diag(sns.histplot, hue=None)
    g.set(rasterized=True)
    g.add_legend()

    # Retrieve the underlying figure from the PairGrid
    fig = g.fig
    return fig


def plot_hist(data, key, y_label):
    """Plot distribution of unique values.
    data - dataframe
    key - key in data, distribution of unique values of which is concidering
    y_label - label of y axe in the plot, describing what the key object is
    """

    hits_number_per_object = defaultdict(int)
    for key_object in data[key]:
        hits_number_per_object[key_object] += 1
    key_objects_cnt_per_hits_number = defaultdict(int)
    for hits_number in hits_number_per_object.values():
        key_objects_cnt_per_hits_number[hits_number] += 1
    N_hits = range(max(key_objects_cnt_per_hits_number))
    N_key_objects = [key_objects_cnt_per_hits_number[i] for i in N_hits]
    fig, ax = plt.subplots()
    ax.set_rasterized(True)
    ax.bar(N_hits, N_key_objects)
    ax.set_xlabel('Number of hits')
    ax.set_ylabel(y_label)
    ax.set_xlim(0, len(N_hits) + 1)

    return fig


def plot_hits_in_taxa(taxon_counts):
    taxon_values = taxon_counts.index.values.tolist()
    taxon_counts = list(taxon_counts)

    fig, ax = plt.subplots(1, 1, figsize=(9, 6), sharey='all')
    ax.bar(taxon_values, taxon_counts)
    ax.set_rasterized(True)
    ax.set_yscale('log')

    ax.set_title('Total hits in taxa')
    ax.set_ylabel('N hits')
    ax.set_xlabel('Taxon')

    for label in ax.get_xticklabels():  # rotate ticks
        label.set_rotation(40)
        label.set_horizontalalignment('right')

    plt.tight_layout()

    return fig


def plot_hits_per_genome_per_taxon(df_handle):
    # plot average N_hits per genome in phyla
    genome_distribution = df_handle['assembly'].value_counts(sort=True).to_dict()
    genome_df = df_handle[['assembly', 'phylum']].drop_duplicates()
    genome_df = genome_df.replace({'assembly': genome_distribution})
    genome_mean = genome_df.groupby('phylum').mean()

    phylum_values = genome_mean.index.values.tolist()
    taxon_mean_hits = genome_mean['assembly'].tolist()

    fig, ax = plt.subplots(1, 1, figsize=(9, 6), sharey=True)
    ax.bar(phylum_values, taxon_mean_hits)
    ax.set_rasterized(True)

    for label in ax.get_xticklabels():  # rotate ticks
        label.set_rotation(40)
        label.set_horizontalalignment('right')

    plt.tight_layout()

    ax.set_title('Average N hits in a genome')
    ax.set_ylabel('N hits')
    ax.set_xlabel('Phylum')

    return fig


def plot_annotation_distribution(df_handle):
    # plot annotation distribution
    annotation_distribution = df_handle['product'].value_counts(sort=True)
    products = annotation_distribution.index.values.tolist()
    product_counts = list(annotation_distribution)
    filter_map = []
    i = 0
    cnt = 0
    number_of_hsps = len(df_handle.index)
    while cnt < number_of_hsps * 0.9:
        filter_map.append(products[i])
        cnt += product_counts[i]
        i += 1

    df_handle['function'] = df_handle['product'].apply(filter_top_column_values, args=[
        filter_map])  # substitute products to "Other" according f_map
    function_distribution = df_handle['function'].value_counts(sort=True)
    functions = function_distribution.index.values.tolist()
    function_distribution = list(function_distribution)

    # put 'Other' to the end of `functions` if present
    put_other_to_the_end(functions, function_distribution)

    fig, ax = plt.subplots(1, 1, figsize=(9, 6), sharey='all')
    ax.bar(functions, function_distribution)
    ax.set_title('Distribution of functional annotations')
    ax.set_xlabel('Product')
    ax.set_ylabel('N hits')
    for label in ax.get_xticklabels():  # rotate ticks
        label.set_rotation(40)
        label.set_horizontalalignment('right')

    return fig


def find_overlaps(starts, ends, length):
    """Return distribution of hits pairwise 
    overlaps lengthes in one genome.
    starts - array with coordinates of hits starts
    ends - array with coordinates of hits ends
    length - array with hits lengthes

    TODO: rewrite to O(n) algorithm 
    """

    hits_pairs_number_per_overlap_length = defaultdict(int)
    n_hits = len(starts)
    genome_length = -1
    for i in range(n_hits):
        flag = False
        start_i = starts[i]
        end_i = ends[i]
        if start_i > end_i:
            if genome_length == -1:
                genome_length = length[i] + start_i - end_i
            end_i += genome_length
            flag = True
        for j in range(n_hits):
            if i == j:
                continue
            start_j = starts[j]
            end_j = ends[j]
            if flag:
                if start_j < end_j:
                    start_j += genome_length
                end_j += genome_length
            if start_i < start_j < end_j < end_i or start_i < start_j - genome_length < end_j - genome_length < end_i:
                hits_pairs_number_per_overlap_length[length[j]] += 1
            elif start_i < start_j < end_i:
                hits_pairs_number_per_overlap_length[end_i - start_j] += 1
    return hits_pairs_number_per_overlap_length


def plot_overlaps(data):
    """Plot distribution of hits pairwise overlaps lengthes.
    data - dataframe
    """

    hits_pairs_number_per_overlap_length = defaultdict(int)
    for assembly in data['assembly'].unique():
        starts = data[data['assembly'] == assembly]['start']
        ends = data[data['assembly'] == assembly]['end']
        length = data[data['assembly'] == assembly]['protein_length']
        for overlap_length, hits_pairs_number in find_overlaps(starts, ends, length).items():
            hits_pairs_number_per_overlap_length[overlap_length] += hits_pairs_number
    overlap_lengths = range(1, max(hits_pairs_number_per_overlap_length) + 1)
    N_hits_pairs = [hits_pairs_number_per_overlap_length[i] for i in overlap_lengths]
    fig, ax = plt.subplots()
    ax.set_rasterized(True)
    ax.bar(overlap_lengths, N_hits_pairs)
    ax.set_xlabel('Overlap length')
    ax.set_ylabel('Number of hits pairs')
    ax.set_xlim(0, len(overlap_lengths) + 1)

    return fig


def plot_3d(df, x_col, y_col, z_col, color_axis, color_dict=None):
    """Plot 3d scatter plot.
    data_points - list of tuples (x, y, z)
    df - dataframe with data for annotation
    color_axis - column that defines colors of the dots
    color_dict
    """

    xs = list(df[x_col])
    ys = list(df[y_col])
    zs = list(df[z_col])

    data_points = [(x, y, z) for x, y, z in zip(xs, ys, zs)]

    if color_dict is not None:
        print(color_dict)
        print(list(df[color_axis]))
        colors = [color_dict[value] for value in list(df[color_axis])]
    else:
        colors = [value for value in list(df[color_axis])]
    labels = [value for value in list(df[color_axis])]
    markers = ['x' if q == 'plasmid' else 'o' for q in list(df['replicon_type'])]

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize=(9, 6))
    ax.set_rasterized(True)

    # construct legend
    if color_dict is not None:
        legend_elements = [Line2D([0], [0], marker='o', color='w',
                                  label=value, markerfacecolor=color_dict[value])
                           for value in df[color_axis].unique().tolist()]
    else:
        legend_elements = None

    # plot data
    if color_dict is not None:
        for point, color, mark, lab in zip(data_points, colors, markers, labels):
            x, y, z = point
            ax.scatter(x, y, z, c=color, marker=mark, label=lab,
                       alpha=0.4, s=3, linewidth=0.8)
    else:
        xs = []
        ys = []
        zs = []
        for x, y, z in data_points:
            xs.append(x)
            ys.append(y)
            zs.append(z)
        sc = ax.scatter(xs, ys, zs, c=colors, norm=matplotlib.colors.PowerNorm(gamma=0.1),
                        edgecolors='black', marker='o', alpha=0.8, s=3, linewidth=0.2, cmap='jet')

    # set upper limits
    ax.set_xlim(None, ident_lim)
    ax.set_ylim(None, overlap_lim)
    ax.set_zlim(None, length_lim)

    # set axis labels
    ax.set_xlabel('Log e-value')
    ax.set_ylabel('Query coverage, %')
    ax.set_zlabel('Hit length, aa')

    if color_dict is not None:
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.1, 1), loc='upper left')
    else:
        plt.colorbar(sc, label='E-value',
                     ticks=[10, 1, 0.1, 0.05, 0.01, 0.001, 0.0001, 10 ** -6, 10 ** -9, 10 ** -12, 10 ** -100],
                     format='%.0e', pad=0.1)
    plt.tight_layout()

    return fig


if __name__ == '__main__':
    try:
        # arguments
        project = sys.argv[1]  # argument -- project name
        database = sys.argv[2]  # arument -- database name

        # upper limits for 3D plots - optional arguments
        try:
            ident_lim = int(sys.argv[3])
        except:
            ident_lim = None

        try:
            overlap_lim = int(sys.argv[4])
        except:
            overlap_lim = None

        try:
            length_lim = int(sys.argv[5])
        except:
            length_lim = None

        # log start to exit log
        exitlog_path = Path('projects') / project / 'exit_log.txt'
        with open(exitlog_path, 'a') as outfile:
            outfile.write('3 started\n')

        # read data
        data = pd.read_csv(Path('projects') / project / 'hits_df.csv',
                           index_col=0)  # dataframe with hmmer results

        # pre-filtering of hits with e-value > 0.1
        data = data[data.evalue <= 0.1]

        data = filter_top_column_values(df=data, filter_column='phylum', assign_column='taxon', keep_first=12)
        data = filter_top_column_values(df=data, filter_column='product', assign_column='function', keep_first=12)

        # output of information about database, query, homology search tool, and  
        # numerical characteristics of search results

        hmmer_results_path = Path('projects') / project / 'hits.txt'
        hmmer_results = SearchIO.read(hmmer_results_path, "hmmer3-text")

        # variables needed to generate plots

        taxon_counts = data['taxon'].value_counts(sort=True)
        taxon_counts = put_other_to_end(taxon_counts)

        print(type(taxon_counts))
        print(taxon_counts)
        print(taxon_counts.keys())
        taxon_colormap = generate_colormap(taxon_counts.index.tolist())

        function_counts = data['function'].value_counts(sort=True)
        function_counts = put_other_to_end(function_counts.keys())

        function_colormap = generate_colormap(function_counts.index.tolist())

        # list of matplotlib plots to generate and export to pdf

        text_report = generate_text_report(data, hmmer_results)
        plots = [
            plot_hist(data, 'assembly', 'Number of genomes'),
            plot_hist(data, 'lcs', 'Number of proteins'),
            # plot_overlaps(data),  # plot hits overlaps lengths distribution
            plot_3d(df=data,
                    x_col='lg_evalue',
                    y_col='query_coverage',
                    z_col='protein_length',
                    color_axis='taxon',
                    color_dict=taxon_colormap
                    ),
            pairplot(data=data,
                     columns=['lg_evalue', 'query_coverage', 'protein_length', 'taxon'],
                     hue='taxon',
                     palette=taxon_colormap
                     ),
            plot_hits_in_taxa(taxon_counts),
            plot_hits_per_genome_per_taxon(data),
            plot_annotation_distribution(data),
            plot_3d(df=data,
                    x_col='lg_evalue',
                    y_col='query_coverage',
                    z_col='protein_length',
                    color_axis='function',
                    color_dict=function_colormap
                    )
        ]

        pdf_path = Path('projects') / project / 'reports' / 'homology_search_report.pdf'
        save_to_pdf(text_report, plots, pdf_path)

    except Exception as e:
        ecx_type = str(type(e))
        print(e)
        print(f"EXCEPTION {ecx_type} RAISED")
        with open(exitlog_path, 'a') as outfile:
            outfile.write('3 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 nto the exit log
        with open(exitlog_path, 'a') as outfile:
            outfile.write('3 0\n')

