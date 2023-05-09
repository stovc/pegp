"""Plot results of hits filtration and clusterization

- input:  "hits_df.csv"
- output: "filtering_clustering_report.pdf"

Output contains folowing:
    Hits filtering:
        filtering parameters - evalue, query coverage
        number of input hits
        number of output hits
    Hits clusterization:
        clusterization programm
        clusterization parameter (currently 90)
        number of input hits
        number of output hits
        number of hits after genome consistency procedure 
        3d plot of log_evalue, query coverage, length colored by:
            - not pass filtering
            - pass filtering but not representative
            - representative
            - return after genome consistency procedure

"""

import sys
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import subprocess
import traceback


def plot_3d(data_points, df, color_axis, color_dict=None):
    """Plot 3d scatter plot.
    data_points - list of tuples (x, y, z)
    df - dataframe with data for annotation
    color_axis - column that defines colors of the dots
    color_dict"""
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
    ax.set_xlabel('Log e-value')
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
        
        exitlog_path = Path('projects') / project / 'exit_log.txt'
        with open(exitlog_path, 'a') as outfile:
            outfile.write('8 started\n')

        # read data
        data_path = Path('projects') / project / 'hits_df.csv'  # path to the dataframe with hits
        data = pd.read_csv(data_path, index_col=0)  # dataframe with hits
        filtered_data = data[data.filtered == True]
        filtered_clustered_data = data[data.filtered_clustered == True]

        # create output to plot to
        out_path = Path('projects') / project / 'filtering_clustering_report.pdf'  # path to the output report pdf
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_path)  # report pdf object

        # output information about filtering results
        page = plt.figure(figsize=(11.69,8.27))
        page.clf()
        txt = "Hits filtering report"
        page.text(0.1, 0.9, txt, transform=page.transFigure, size=24, ha="left")

        # filtering parameters
        txt = "Filtering parameters:\n"
        filtering_log_path = Path('projects') / project / 'filtering_log.txt'
        with open(filtering_log_path, 'r') as filtering_log:
            for line in filtering_log:
                txt += ' - ' + line + '\n'
        page.text(0.1, 0.55, txt, transform=page.transFigure, size=20, ha="left")

        # number of input hits
        filtering_input_hits_cnt = len(data.index)
        txt = f"Number of filtering input hits: {filtering_input_hits_cnt}"
        page.text(0.1, 0.5, txt, transform=page.transFigure, size=20, ha="left")

        # number of output hits
        filtering_output_hits_cnt = len(filtered_data.index)
        txt = f"Number of filtering output hits: {filtering_output_hits_cnt}"
        page.text(0.1, 0.4, txt, transform=page.transFigure, size=20, ha="left")

        pdf.savefig()

        # output information about clustering and genome consistancy 
        page = plt.figure(figsize=(11.69,8.27))
        page.clf()
        txt = "Hits clustering and genome consistancy report"
        page.text(0.1, 0.9, txt, transform=page.transFigure, size=24, ha="left")

        # clustering programm
        clustering_log_path = Path('projects') / project / 'clustering_log.txt'
        with open(clustering_log_path, 'r') as clustering_log:
            version_line = [line for line in clustering_log][0]
            version_list = version_line.split()
            version = ' '.join(version_list[1:-1])
        txt = f"Clustering program version: {version}"
        page.text(0.1, 0.8, txt, transform=page.transFigure, size=20, ha="left") 

        # clustering parameter
        clustering_parameter = 90 # this is hardcoded now!
        txt = f"Clustering parameter: {clustering_parameter}"
        page.text(0.1, 0.7, txt, transform=page.transFigure, size=20, ha="left")       

        # number of input hits
        clustering_input_hits_cnt = filtering_output_hits_cnt
        txt = f"Number of clustering input hits: {clustering_input_hits_cnt}"
        page.text(0.1, 0.6, txt, transform=page.transFigure, size=20, ha="left")

        # number of output hits
        clustered_faa_path = Path('projects') / project / 'clustered90.faa'
        with open(clustered_faa_path) as file:
            clustering_output_hits_cnt = sum([1 for record in SeqIO.parse(file, 'fasta')])
        txt = f"Number of clustering output hits: {clustering_output_hits_cnt}"
        page.text(0.1, 0.5, txt, transform=page.transFigure, size=20, ha="left")

        # number of hits after genome consistency
        genome_consistency_output_hits_cnt = len(filtered_clustered_data.index)
        txt = f"Number of hits after genome consistency procedure: {genome_consistency_output_hits_cnt}"
        page.text(0.1, 0.4, txt, transform=page.transFigure, size=20, ha="left") 

        pdf.savefig()

        #3d plot

        pdf.close()


    except Exception as e:
        ecx_type = str(type(e))
        print(f"EXCEPTION {ecx_type} RAISED")

        with open(exitlog_path, 'a') as outfile:
            outfile.write('8 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        with open(exitlog_path, 'a') as outfile:
            outfile.write('8 0\n')
