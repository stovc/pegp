"""Plot results of triming

- input:  "aligned.fa", "trimed.fa"
- output: "triming_report.pdf"

Output contains folowing information:
    - triming program
    - triming threshold
    - number of columns before triming
    - number of columns after triming

"""

import sys  
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import traceback


if __name__ == '__main__':
    try:
        # arguments
        project = sys.argv[1]  # argument -- project name
        
        exitlog_path = Path('projects') / project / 'exit_log.txt'
        with open(exitlog_path, 'a') as outfile:
            outfile.write('8 started\n')

        # read data
        aligned_path = Path('projects') / project / 'aligned.fa'  # path to fasta file with aligned reads before trimming
        trimed_path = Path('projects') / project / 'trimed.fa'  # path to fasta file with aligned reads after trimming

        # create output to plot to
        out_path = Path('projects') / project / 'triming_report.pdf'  # path to the output report pdf
        pdf = matplotlib.backends.backend_pdf.PdfPages(out_path)  # report pdf object

        # create a page to output information
        page = plt.figure(figsize=(11.69,8.27))
        page.clf()
        txt = "Triming report"
        page.text(0.1, 0.9, txt, transform=page.transFigure, size=24, ha="left")

        # get information from triming_log.txt
        triming_log_path = Path('projects') / project / 'triming_log.txt'
        with open(triming_log_path, 'r') as f:
            triming_log = [line for line in f]

        # trimming program
        txt = f"Triming program: {triming_log[1]}"
        page.text(0.1, 0.8, txt, transform=page.transFigure, size=20, ha="left")

        # triming threshold 
        txt = f"Triming threshold: {triming_log[3]}"
        page.text(0.1, 0.7, txt, transform=page.transFigure, size=20, ha="left")

        # number of columns before triming
        aligned_seq = [record.seq for record in SeqIO.parse(aligned_path, 'fasta')]
        aligned_len = len(aligned_seq[0])
        txt = f"Number of columns before triming: {aligned_len}"
        page.text(0.1, 0.65, txt, transform=page.transFigure, size=20, ha="left")

        # number of columns after triming
        trimed_seq = [record.seq for record in SeqIO.parse(trimed_path, 'fasta')]
        trimed_len = len(trimed_seq[0])
        txt = f"Number of columns after triming: {trimed_len}"
        page.text(0.1, 0.55, txt, transform=page.transFigure, size=20, ha="left")

        pdf.savefig()

        pdf.close()

    except Exception as e:
        ecx_type = str(type(e))
        print(e)
        print(f"EXCEPTION {ecx_type} RAISED")

        with open(exitlog_path, 'a') as outfile:
            outfile.write('8 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        with open(exitlog_path, 'a') as outfile:
            outfile.write('8 0\n')
