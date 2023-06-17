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
"""

import sys
from pathlib import Path
from Bio import SeqIO, SearchIO
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.backends.backend_pdf
import seaborn as sns
import subprocess
import traceback


def filter_df(object, filter_map):
    """Return object if it belongs to the filter map. Otherwise, return 'other'. 
    """

    if object in filter_map:
        return object
    else:
        return 'Other'


def put_other_to_the_end(data, distribution):
    """Put 'Other' to the end of `data` if present
    data - array with features of hits
    distribution - number of hits with particular feature
    """
    
    if 'Other' in data:
        other_position = data.index('Other')
        other_count = distribution[other_position]
        del data[other_position]
        del distribution[other_position]
        data.append('Other')
        distribution.append(other_count)


def pairplot(data, columns, hue, palette):
    """Plot seaborn pairplot.
    Data - dataframe
    Columns - columns from data to be plotted
    Hue - column that defines colors of plotted dots
    Palette - color palette fot dots
    """
    g = sns.PairGrid(data[columns], diag_sharey=False, hue=hue,
                     height=3, aspect=1.5, palette=palette, corner=True)
    g.map_lower(sns.scatterplot, edgecolor="k", linewidth=0.2, s=3)
    g.map_diag(sns.histplot, hue=None)
    g.set(rasterized=True)

    pdf.savefig(dpi=300, bbox_inches='tight')

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
    ax.set_xlim(0, len(N_hits)+1)
    pdf.savefig(dpi=300, bbox_inches='tight')


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
        length = data[data['assembly'] == assembly]['length']
        for overlap_length, hits_pairs_number in find_overlaps(starts, ends, length).items():
            hits_pairs_number_per_overlap_length[overlap_length] += hits_pairs_number
    overlap_lengths = range(1, max(hits_pairs_number_per_overlap_length) + 1)
    N_hits_pairs = [hits_pairs_number_per_overlap_length[i] for i in overlap_lengths]
    fig, ax = plt.subplots()
    ax.set_rasterized(True)
    ax.bar(overlap_lengths, N_hits_pairs)
    ax.set_xlabel('Overlap length')
    ax.set_ylabel('Number of hits pairs')
    ax.set_xlim(0, len(overlap_lengths)+1)
    pdf.savefig(dpi=300, bbox_inches='tight')    


def plot_3d(data_points, df, color_axis, color_dict=None):
    """Plot 3d scatter plot.
    data_points - list of tuples (x, y, z)
    df - dataframe with data for annotation
    color_axis - column that defines colors of the dots
    color_dict
    """

    if color_dict is not None:
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
        plt.colorbar(sc, label='E-value', ticks=[10, 1, 0.1, 0.05, 0.01, 0.001, 0.0001, 10**-6, 10**-9, 10**-12, 10**-100],
                     format='%.0e', pad=0.1)
    plt.tight_layout()

    pdf.savefig(dpi=400, bbox_inches='tight')


if __name__ == '__main__':
    try:
        # arguments
        project = sys.argv[1]  # argument -- project name
        database = sys.argv[2] # arument -- database name

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
        data_path = Path('projects') / project / 'hits_df.csv'  # path to the hmmer result dataframe
        data = pd.read_csv(data_path, index_col=0)  # dataframe with hmmer results
       
        # create output to plot to
        out_path = Path('projects') / project / 'search_report.pdf'  # path to the output report pdf
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_path)  # report pdf object
        
        # set of colors for plots
        COLORS20 = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4',
                    '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8',
                    '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9'] + 30 * ['#a9a9a9']
        
        #pre-filtering of hits with e-value > 0.1
        df_handle = data[data.evalue <= 0.1]
        
        df_handle = df_handle.replace({'replicon_type': {'main': 'chromosome'}})

        # output of information about database, query, homology search tool, and  
        # numerical characteristics of search results
        firstPage = plt.figure(figsize=(11.69,8.27))
        firstPage.clf()

        txt = "Preliminary filtering was carried out by e-value < 0.1"
        firstPage.text(0.1, 0.9, txt, transform=firstPage.transFigure, size=20, ha="left")

        txt = f"Database: {database}"
        firstPage.text(0.1, 0.8, txt, transform=firstPage.transFigure, size=20, ha="left")

        protein_faa_path = Path('databases') / database / 'protein.faa'
        with open(protein_faa_path) as file:
            database_size = sum([1 for record in SeqIO.parse(file, 'fasta')])
        txt = f"Database size: {str(database_size)}"
        firstPage.text(0.1, 0.7, txt, transform=firstPage.transFigure, size=20, ha="left")

        hmmer_results_path = Path('projects') / project / 'hits.txt'
        hmmer_results = SearchIO.read(hmmer_results_path, "hmmer3-text")
        txt = f"Query: {hmmer_results._id}"
        firstPage.text(0.1, 0.6, txt, transform=firstPage.transFigure, size=20, ha="left")

        txt = f"Homology search using: {hmmer_results.program} {hmmer_results.version}"
        firstPage.text(0.1, 0.5, txt, transform=firstPage.transFigure, size=20, ha="left")

        number_of_HSPs = len(df_handle.index)
        txt = f"Number of HSPs: {str(number_of_HSPs)}"
        firstPage.text(0.1, 0.4, txt, transform=firstPage.transFigure, size=20, ha="left")

        hit_proteins_number = len(df_handle['lcs'].unique())
        txt = f"Number of hit proteins: {str(hit_proteins_number)}"
        firstPage.text(0.1, 0.3, txt, transform=firstPage.transFigure, size=20, ha="left")   

        hit_genomes_number = len(df_handle['assembly'].unique())
        txt = f"Number of hit genomes: {str(hit_genomes_number)}"
        firstPage.text(0.1, 0.2, txt, transform=firstPage.transFigure, size=20, ha="left")      

        pdf.savefig()

        # plot genomes distribution
        plot_hist(df_handle, 'assembly', 'Number of genomes')
        
        # plot proteins distribution
        plot_hist(df_handle, 'lcs', 'Number of proteins')
       
        # plot hits overlaps lengths distribution
        # plot_overlaps(df_handle)

        # 3d plot of identity, query_coverage, and length colored by taxon
        # Taxon distribution
        # Finally, generates Taxons list and Taxon_counts list
        taxon_destribution = df_handle['phylum'].value_counts(sort=True)

        taxons = taxon_destribution.index.values.tolist()

        filter_map = taxons[0:13]  # selects 1st N taxons
        df_handle['taxon'] = df_handle['phylum'].apply(filter_df, args=[filter_map])  # substitute taxons to "other" according f_map
        taxon_destribution = df_handle['taxon'].value_counts(sort=True)

        taxons = taxon_destribution.index.values.tolist()
        taxon_counts = list(taxon_destribution)

        # put 'Other' to the end of `taxons` if present
        put_other_to_the_end(taxons, taxon_counts)

        # generate taxon color dicts
        taxon_color_dict = {}
        for i in range(len(taxons)):
            taxon_color_dict[taxons[i]] = COLORS20[i]

        # 3d plot
        xs = list(df_handle['lg_evalue'])
        ys = list(df_handle['query_coverage'])
        zs = list(df_handle['protein_length'])

        data_points = [(x, y, z) for x, y, z in zip(xs, ys, zs)]
        plot_3d(data_points, df_handle, color_axis='taxon', color_dict=taxon_color_dict)

        # # 3d plot of identity, query_coverage, and length colored by evalue
        # plot_3d(data_points, df_handle, color_axis='lg_evalue')

        # plot pairplot with taxonomy information
        cols = ['lg_evalue', 'query_coverage', 'protein_length', 'taxon']  # 'evalue^0.1', 'length', 'identity',
        pairplot(df_handle, cols, 'taxon', taxon_color_dict)

      
        # plot taxon distribution
        fig, ax = plt.subplots(1, 1, figsize=(9, 6), sharey='all')
        ax.bar(taxons, taxon_counts)
        ax.set_rasterized(True)
        ax.set_yscale('log')

        ax.set_title('Total hits in phyla')
        ax.set_ylabel('N hits')
        ax.set_xlabel('Phylum')

        for label in ax.get_xticklabels():  # rotate ticks
            label.set_rotation(40)
            label.set_horizontalalignment('right')

        plt.tight_layout()
        pdf.savefig(bbox_inches='tight')
       
        # plot average N_hits per genome in phyla
        genome_distribution = df_handle['assembly'].value_counts(sort=True).to_dict()
        genome_df = df_handle[['assembly', 'phylum']].drop_duplicates()
        genome_df = genome_df.replace({'assembly': genome_distribution})
        genome_mean = genome_df.groupby('phylum').mean()

        taxons = genome_mean.index.values.tolist()
        taxon_mean_hits = genome_mean['assembly'].tolist()

        fig, ax = plt.subplots(1, 1, figsize=(9, 6), sharey=True)
        ax.bar(taxons, taxon_mean_hits)
        ax.set_rasterized(True)

        for label in ax.get_xticklabels():  # rotate ticks
            label.set_rotation(40)
            label.set_horizontalalignment('right')

        plt.tight_layout()

        ax.set_title('Average N hits in a genome')
        ax.set_ylabel('N hits')
        ax.set_xlabel('Phylum')

        pdf.savefig(bbox_inches='tight')

       # plot annotation distribution
        annotation_destribution = df_handle['product'].value_counts(sort=True)
        products = annotation_destribution.index.values.tolist()
        product_counts = list(annotation_destribution)
        filter_map = []
        i = 0
        cnt = 0
        while cnt < number_of_HSPs * 0.9:
            filter_map.append(products[i])
            cnt += product_counts[i]
            i += 1
        
        df_handle['function'] = df_handle['product'].apply(filter_df, args=[filter_map])  # substitute products to "Other" according f_map
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
        pdf.savefig(bbox_inches='tight')
        
        # output the list of all functional annotations
        i = 0
        while i < len(products):
            page = plt.figure(figsize=(11.69,8.27))
            page.clf()
            if i == 0:
                page.text(0.1, 0.9, "List of all functional annotations", transform=page.transFigure, size=24, ha="left")
            txt = ""
            until = i + 25
            while i < len(products) and i < until:
                txt += products[i] + ':   ' + str(product_counts[i]) + '\n'
                i += 1
            page.text(0.1, 0.1, txt, transform=page.transFigure, size=15, ha="left")
            pdf.savefig()
        
        # 3d plot of identity, query_coverage, and length colored by functional annotation
        function_color_dict = {}
        for i in range(len(functions)):
            function_color_dict[functions[i]] = COLORS20[i]
        plot_3d(data_points, df_handle, color_axis='function', 
                                        color_dict=function_color_dict)

        pdf.close()
        
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
