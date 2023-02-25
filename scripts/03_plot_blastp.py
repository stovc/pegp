"""Plot BLAST hits and their properties.

- Step 3 in the pipeline
"""

import sys
import getopt
from pathlib import Path
import numpy as np
from math import log2
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors, cm
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
import matplotlib.backends.backend_pdf
import seaborn as sns
import subprocess
import traceback


def filer_taxonomy(taxon, filter_map):
    """Return taxon if it belongs to the filter map. Otherwise, return 'other'. """
    if taxon in filter_map:
        return taxon
    else:
        return 'Other'


def pairplot(data, columns, hue, palette):
    g = sns.PairGrid(data[columns], diag_sharey=False, hue=hue,
                     height=3, aspect=1.5, palette=palette, corner=True)
    g.map_lower(sns.scatterplot, edgecolor="k", linewidth=0.2, s=3)
    g.map_diag(sns.histplot, hue=None)
    g.set(rasterized=True)

    pdf.savefig(dpi=300)


def plot_3d(data_points, df, color_axis, color_dict):
    """Plot 3d scatter plot."""
    if color_dict is not None:
        colors = [color_dict[value] for value in list(df[color_axis])]
    else:
        colors = [value for value in list(df[color_axis])]
    labels = [value for value in list(df_handle[color_axis])]
    markers = ['x' if q == 'plasmid' else 'o' for q in list(df['replicon_type'])]

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
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
    ax.set_xlabel('Identity')
    ax.set_ylabel('Query coverage, %')
    ax.set_zlabel('Hit length, aa')

    if color_dict is not None:
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        plt.colorbar(sc, label='E-value', ticks=[10, 1, 0.1, 0.05, 0.01, 0.001, 0.0001, 10**-6, 10**-9, 10**-12, 10**-100],
                     format='%.0e')
    plt.tight_layout()

    pdf.savefig(dpi=400)


if __name__ == '__main__':
    try:
        # arguments
        project = sys.argv[1]  # argument -- project name

        # upper limits for 3D plots - optional arguments
        try:
            ident_lim = int(sys.argv[2])
        except:
            ident_lim = None

        try:
            overlap_lim = int(sys.argv[3])
        except:
            overlap_lim = None

        try:
            length_lim = int(sys.argv[4])
        except:
            length_lim = None

        # generate paths
        data_path = Path('projects') / project / 'blastp_df.csv'  # path to the blastp result dataframe
        out_path = Path('projects') / project / 'blastp_report.pdf'  # path to the output report pdf
        exitlog_path = Path('projects') / project / 'exit_log.txt'

        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '3 started'], stdout=outfile)

        data = pd.read_csv(data_path, index_col=0)  # dataframe with blastp results
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_path)  # report pdf object

        colors20 = ['#e6194B', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#42d4f4',
                    '#f032e6', '#bfef45', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#fffac8',
                    '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9'] + 30*['#a9a9a9']

        # Taxon distribution
        # Finally, generates Taxons list and Taxon_counts list
        df_handle = data
        taxon_destribution = df_handle['phylum'].value_counts(sort=True)

        taxons = taxon_destribution.index.values.tolist()

        filter_map = taxons[0:13]  # selects 1st N taxons
        df_handle['taxon'] = df_handle['phylum'].apply(filer_taxonomy, args=[filter_map])  # substitute taxons to "other" according f_map
        taxon_destribution = df_handle['taxon'].value_counts(sort=True)

        taxons = taxon_destribution.index.values.tolist()
        taxon_counts = list(taxon_destribution)

        other_position = taxons.index('Other')
        other_count = taxon_counts[other_position]

        del taxons[other_position]
        del taxon_counts[other_position]

        taxons.append('Other')
        taxon_counts.append(other_count)

        df_handle['evalue^0.1'] = df_handle['evalue']**0.1
        df_handle = df_handle.replace({'replicon_type': {'main': 'chromosome'}})

        # generate taxon color dicts

        taxon_color_dict = {}
        for i in range(len(taxons)):
            taxon_color_dict[taxons[i]] = colors20[i]

        color_dict_cluster = {'0': 'red', '1': 'blue', '2': 'green', '3': 'orange', '4': 'magenta', '5': 'cyan',
                              '-': 'grey'}


        # 3d plot

        xs = list(df_handle['identity'])
        ys = list(df_handle['query_coverage'])
        zs = list(df_handle['length'])

        data_points = [(x, y, z) for x, y, z in zip(xs, ys, zs)]

        plot_3d(data_points, df_handle, 'taxon', taxon_color_dict)  # colored by taxon

        # 3d plot clusters

        # plot_3d(data_points, df_handle, 'cluster', color_dict_cluster)  # colored by cluster

        # 3d plot evalue
        plot_3d(data_points, df_handle, 'evalue', None)  # colored by evalue

        # plot pairplot with taxonomy information
        cols = ['query_coverage', 'identity', 'length', 'evalue^0.1', 'taxon']
        pairplot(df_handle, cols, 'taxon', taxon_color_dict)

        # plot taxon distribution

        fig, ax = plt.subplots(1, 1, figsize=(9, 6), sharey='all')
        ax.bar(taxons, taxon_counts)
        ax.set_rasterized(True)
        ax.set_yscale('log')

        ax.set_title('Total hits in phyla')
        ax.set_xlabel('N hits')
        ax.set_ylabel('Phylum')

        for label in ax.get_xticklabels():  # rotate ticks
            label.set_rotation(40)
            label.set_horizontalalignment('right')

        plt.tight_layout()
        pdf.savefig()

        # plot average N_hits per genome in phyla
        genome_destribution = df_handle['assembly'].value_counts(sort=True).to_dict()
        genome_df = df_handle[['assembly', 'phylum']].drop_duplicates()
        genome_df = genome_df.replace({'assembly': genome_destribution})
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
        ax.set_xlabel('N hits')
        ax.set_ylabel('Phylum')

        pdf.savefig()

        ###

        pdf.close()

    except Exception as e:
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '3 ' + ecx_type], stdout=outfile)

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        # print exin code 0 nto the exit log
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '3 0'], stdout=outfile)
