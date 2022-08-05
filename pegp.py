import os
import subprocess
from pathlib import Path
import pandas as pd
from prettytable import PrettyTable
from pathlib import Path
from colorama import Fore, Back, Style


# INPUT OUTPUT
def highlight(s):
    s = Fore.BLACK + Back.WHITE + s + Style.RESET_ALL
    return s


def update_status():
    data = pd.DataFrame(DATA, index=list(range(1, 16)))

    for project in projects:

        status_path = PROJECTS_PATH / project / 'status.txt'
        status = pd.read_csv(status_path, sep='\t', index_col=0)

        try:
            exit_log = open(PROJECTS_PATH / project / 'exit_log.txt', 'r')
            new_exit_log = ''

            for line in exit_log:
                step = int(line.split(' ')[0])              # first word
                step_status = line.split(' ')[1][:-1]  # second word without end line symbol \n
                if step_status == 'started':
                    status.at[step, project] = 'running'
                elif step_status.isnumeric():
                    if step_status == '0':
                        status.at[step, project] = 'done'
                        for i in UNLOCKS[step]:
                            if status.at[i, project] == '-':
                                status.at[i, project] = 'ready'
                        print('should make it ready!')

                    else:
                        status.at[step, project] = f'failed {step_status}'
                else:
                    new_exit_log += line
            exit_log.close()

            exit_log = open(PROJECTS_PATH / project / 'exit_log.txt', 'w')
            exit_log.write(new_exit_log)
            exit_log.close()

            status.to_csv(status_path, sep='\t')
        except:
            #  print('exit log open exception')
            ''
        data = pd.concat([data, status], axis=1)
    return data


def draw_table(data):
    table = PrettyTable(['Step'] + projects)
    table.align["Step"] = "l"

    for row in data.iterrows():
        lst_row = list(row[1])
        updated_row = []
        for j in lst_row:
            nj = j.replace('done', Fore.GREEN + 'Done' + Style.RESET_ALL)
            nj = nj.replace('ready', Fore.YELLOW + 'Ready' + Style.RESET_ALL)
            nj = nj.replace('queued', Fore.BLUE + Style.BRIGHT + 'Queued' + Style.RESET_ALL)
            nj = nj.replace('running', Fore.CYAN + 'Running' + Style.RESET_ALL)
            nj = nj.replace('failed', Fore.RED + 'Failed' + Style.RESET_ALL)
            updated_row.append(nj)
        table.add_row(updated_row)
    return table


def update_screen():
    data = update_status()
    table = draw_table(data)

    print('Database:', database)
    print(table)
    print(PROMPT)
    #for i in log.history:
    #    print(i)


class Log:
    history = []

    def inp(self, s):
        self.history.append(f'>>> {s}')

    def out(self, s):
        print(s)
        self.history.append(s)


# COMMANDS

# Start process
def start_process(step, project, database, *args):
    script = SCRIPTS[step]
    path = BATCH_SCRIPTS_PATH / script
    status_path = PROJECTS_PATH / project / 'status.txt'
    process = subprocess.run(["sbatch", path, project, database] + list(args), capture_output=True)

    status = pd.read_csv(status_path, sep='\t', index_col=0)
    status.at[step, project] = 'queued'
    status.to_csv(status_path, sep='\t')

    update_screen()
    return process


# Call different actions
def call(action):
    """Decorator. Take input as a string, parse arguments, check..  WIP.."""

    def wrapper(input):
        n_args = action.__code__.co_argcount
        args = tuple(input.split(' ')[1:])
        if len(args) < n_args:              # TODO: it should check if function fits a preset argument range
            log.out(f'{n_args} arguments expected, {len(args)} passed')
        else:
            action(*args)
    return wrapper


@call
def start_all_steps():
    list_of_ready = []
    for col in projects:
        numbers = data.loc[data[col] == 'ready', col].index.tolist()
        for number in numbers:
            list_of_ready.append((number + 1, col))

    for step, project in list_of_ready:
        log.out(f'Starting Step {step} on {project}')
        start_process(step, project, database)
        log.out('Started')


@call
def start_step(project, step, *args):
    log.out(f'Starting Step {step} on {project}')
    start_process(int(step), project, database, *args)
    log.out('Started')


@call
def create_project(project, inp_path):

    new_project_path = PROJECTS_PATH / project
    new_inp_path = new_project_path / 'input.faa'
    if not os.path.exists(new_project_path):
        os.mkdir(new_project_path)
        os.system(f'cp {inp_path} {new_inp_path}')
        status = open(new_project_path / 'status.txt', 'w')
        status.write(STATUS_FILE.format(project))
        status.close()

        projects.append(project)

        update_screen()

        log.out(f'Created a new project {project}')
    else:
        log.out(f'Project with name {project} already exists')


@call
def delete_project(project):
    delete_project_path = PROJECTS_PATH / project
    if os.path.exists(delete_project_path):
        os.system(f'rm -R {delete_project_path}')

        projects.remove(project)

        update_screen()

        log.out(f'Deleted the project {project}')
    else:
        log.out(f'Project with name {project} does not exist')


@call
def update():
    update_screen()


@call
def quitx():
    quit()


# constants
PROJECTS_PATH = Path("projects/")
BATCH_SCRIPTS_PATH = Path("batch_scripts/")
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

PROMPT = highlight('q') + \
         ': quit, ' + \
         highlight('s [P] [S]') + Style.RESET_ALL + \
         ': start Project [P] Step [S], ' + \
         highlight('a') + \
         ': start all ready processes, ' + \
         highlight('n [P] [Pa]') + \
         ': start new project [P] path [Pa], ' + \
         highlight('d [p]') + \
         ': Delete project [p]'

STATUS_FILE = '''\t{}
1\tready
2\t-
3\t-
4\t-
5\t-
6\t-
7\t-
8\t-
9\t-
10\t-
11\t-
12\t-
13\t-
14\t-
15\t-
'''

DATA = {'Step': ['1. Blast',
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
        }

UNLOCKS = {
    1: [2],
    2: [3],
    3: [4],
    4: [5],
    5: [6, 7],
    6: [8, 9, 10],
    8: [11],
    11: [13],
    13: [15],
    9: [12],
    12: [14],
}

# init
database = 'base1'
projects = os.listdir(PROJECTS_PATH)

command = ''
log = Log()

data = update_status()

table = draw_table(data)
update_screen()


class Command:
    def __init__(self, text, action):
        self.text = text
        self.action = action
char = 1

"""
commands = [
    Command(
        text='a',
        action=start_all_steps
    ),
    Command(
        text='s',
        action=start_step
    ),
    Command(
        text='n',
        action=create_project
    ),
    Command(
        text='d',
        action=delete_project
    ),
    Command(
        text='u',
        action=update
    ),
    Command(
        text='q',
        action=quitx
    )
] 
"""

commands = {'a': start_all_steps,
            's': start_step,
            'n': create_project,
            'd': delete_project,
            'u': update,
            'q': quitx}

# input reading
while True:
    command = input('>>> ')
    log.inp(command)

    for i in commands.keys():
        command_args = command.split(' ')
        if command_args[0] == i:
            action = commands[command_args[0]]
            action(command)
            break
    else:
        log.out(f'unknown command "{command}"')
