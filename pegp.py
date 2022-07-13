import os
import subprocess
from pathlib import Path
import pandas as pd
from prettytable import PrettyTable
from pathlib import Path
from colorama import Fore, Back, Style

PROJECTS_LOCATION = Path("projects/")
SLURM_SCRIPTS_LOCATION = Path("slurm_commands/")
SCRIPTS = {1: '01_blastp.sh',
           2: '02_mk_blastp_df.sh',
           3: '03_plot_blastp.sh',
           4: '04_filter_blastp.sh',
           5: '05_cdhit.sh',
           6: '06_mk_gen_consistent.sh',
           7: '07_cluster_dict.sh',
           8: '08_align.sh',
           9: '09_get_translations.sh',
           10: '10_get_gen_context.sh',
           11: '11_trim.sh',
           12: '12_hmmscan.sh',
           13: '13_prottest.sh',
           14: '14_get_domains.sh',
           15: '15_raxml.sh'
           }


def start_process(step, project):
    script = SCRIPTS[step]
    path = SLURM_SCRIPTS_LOCATION / script
    process = subprocess.run(["bash", path, project], shell=True)
    return process


PROJECTS_LOCATION = Path("projects/")

# Prettytable color
R = "\033[0;31;40m" #RED
G = "\033[0;32;40m" # GREEN
Y = "\033[0;33;40m" # Yellow
B = "\033[0;34;40m" # Blue
N = "\033[0m" # Reset


projects = os.listdir(PROJECTS_LOCATION)
data = pd.DataFrame({'Step': ['1. Blast',
                              '2. Blast data',
                              '3. Blast report',
                              '4. Filter',
                              '5. Cluster',
                              '6. Genome consistence',
                              '7. Cluster data',
                              '8. Align',
                              '9. Full translations',
                              '10. Genome context',
                              '11. Trim',
                              '12. Domain search',
                              '13. Prottest',
                              '14. Domain data',
                              '15. Tree']
                     })

for project in projects:
    status = pd.read_csv(PROJECTS_LOCATION / project / 'status.txt', sep='\t',
                         usecols=[2],
                         names=[project]
                         )

    data = pd.concat([data,status], axis=1)

t = PrettyTable(['Step'] + projects)
t.align["Step"] = "l"
for row in data.iterrows():
    lst_row = list(row[1])
    updated_row = []
    for j in lst_row:
        nj = j.replace('done', Fore.GREEN + 'Done' + Style.RESET_ALL)
        nj = nj.replace('ready', Fore.YELLOW + 'Ready' + Style.RESET_ALL)
        nj = nj.replace('running', Fore.CYAN + 'Running' + Style.RESET_ALL)
        updated_row.append(nj)
    t.add_row(updated_row)

print(t)

command = ''
while command != 'q':
    command = input('q: quit, a: start all ready processes. Command: ')

    if command == 'a':
        list_of_ready = []
        for col in projects:
            numbers = data.loc[data[col] == 'ready', col].index.tolist()
            for number in numbers:
                list_of_ready.append((number + 1, col))
        print(list_of_ready)

        for step, project in list_of_ready:
            print(f'Starting Step {step} on {project}')
            process = start_process(step, project)
            print('process', process)
            print('Started')
