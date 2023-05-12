"""Plot results of triming

- input:  "aligned.fa", "trimed.fa", "trim.html"
- output: "triming_report.pdf"

Output contains folowing:
    - triming program
    - triming threshold
    - number of columns before triming
    - number of columns after triming
    - trimAl program report

"""

import sys  
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import traceback
from weasyprint import HTML, CSS
import PyPDF2
import os


if __name__ == '__main__':
    try:
        # arguments
        project = sys.argv[1]  # argument -- project name
        
        exitlog_path = Path('projects') / project / 'exit_log.txt'
        with open(exitlog_path, 'a') as outfile:
            outfile.write('12 started\n')

        # read data
        aligned_path = Path('projects') / project / 'aligned.fa'  # path to fasta file with aligned reads before trimming
        trimed_path = Path('projects') / project / 'trimed.fa'  # path to fasta file with aligned reads after trimming

        # path to the output report pdf
        out_path = Path('projects') / project / 'triming_report.pdf'  

        #create path to temporary pdf file with triming information 
        temp_out_path = Path('projects') / project / 'temp_triming_info_report.pdf'
        pdf = matplotlib.backends.backend_pdf.PdfPages(temp_out_path)  # report pdf object

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
        
        # convert trimAl html report to .pdf
        trimAl_report_path = Path('projects') / project / 'trim.html' 
        trimAl_pdf_report_path = Path('projects') / project / 'trim_temp.pdf' 
        css = CSS(string='''@page {size: 15in 8.27in; margin: 0in 0in 0in 0in;}''')
        trimAl_pdf_report = HTML(trimAl_report_path).write_pdf(stylesheets=[css])
        with open(trimAl_pdf_report_path, 'wb') as f:
            f.write(trimAl_pdf_report)

        # merge triming information file with trimAl report
        merger = PyPDF2.PdfFileWriter()
        temp_out = open(temp_out_path, 'rb')
        temp_out_reader = PyPDF2.PdfFileReader(temp_out)
        for page_num in range(temp_out_reader.numPages):
            page = temp_out_reader.getPage(page_num)
            merger.addPage(page)

        trimAl_pdf_report = open(trimAl_pdf_report_path, 'rb')
        trimAl_pdf_report_reader = PyPDF2.PdfFileReader(trimAl_pdf_report)
        for page_num in range(trimAl_pdf_report_reader.numPages):
            page = trimAl_pdf_report_reader.getPage(page_num)
            merger.addPage(page)

        with open(out_path, 'wb') as f:
            merger.write(f)

        temp_out.close()
        trimAl_pdf_report.close()

        os.remove(temp_out_path)
        os.remove(trimAl_pdf_report_path)

    except Exception as e:
        ecx_type = str(type(e))
        print(e)
        print(f"EXCEPTION {ecx_type} RAISED")

        with open(exitlog_path, 'a') as outfile:
            outfile.write('12 ' + ecx_type + '\n')

        with open('log.txt', 'a') as outfile:
            traceback.print_exc(file=outfile)

    else:
        with open(exitlog_path, 'a') as outfile:
            outfile.write('12 0\n')