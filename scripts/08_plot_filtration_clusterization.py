"""Plot results of hits filtration and clusterization
"""

import sys
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import subprocess
import traceback


if __name__ == '__main__':
    try:
        # arguments
        project = sys.argv[1]  # argument -- project name
        
        exitlog_path = Path('projects') / project / 'exit_log.txt'
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '8 started'], stdout=outfile)

        # read data
        data_path = Path('projects') / project / 'hits_df.csv'  # path to the dataframe with hits
        data = pd.read_csv(data_path, index_col=0)  # dataframe with hits

        # create output to plot to
        out_path = Path('projects') / project / 'filtering_clustering_report.pdf'  # path to the output report pdf
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_path)  # report pdf object

    except Exception as e:
        print("EXCEPTION RAISED")
        ecx_type = str(type(e))

        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '8 ' + ecx_type], stdout=outfile)

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        with open(exitlog_path, 'a') as outfile:
            subprocess.run(["echo", '8 0'], stdout=outfile)
