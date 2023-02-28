import os
import subprocess
from pathlib import Path
import pandas as pd
from prettytable import PrettyTable
from pathlib import Path
from colorama import Fore, Back, Style


# INPUT OUTPUT
def highlight(s):
    """Get string. Return highlighted string."""
    s = Fore.BLACK + Back.WHITE + s + Style.RESET_ALL
    return s


def update_status():
    """Update the status of each step of each project.
    Iterate project directories, look at `status.txt` files.
    Return dataframe with status of steps in projects."""

    data = pd.DataFrame(STEPS, index=list(range(1, 16)))

    for project in projects:

        status_path = PROJECTS_PATH / project / 'status.txt'
        status = pd.read_csv(status_path, sep='\t', index_col=0)

        try:
            exit_log = open(PROJECTS_PATH / project / 'exit_log.txt', 'r')

            for line in exit_log:
                step = int(line.split(' ')[0])  # first word - the step number
                step_status = ''.join(line.split(' ')[1:])[:-1]  # the line except the 1st word and the end line symbol \n - step status
                if step_status == 'started':
                    status.at[step, project] = 'running'

                elif step_status == '0':  # exit code 0 - step successfully completed
                    status.at[step, project] = 'done'

                    # unlock steps dependent on the completed step
                    for i in UNLOCKS[step]:
                        if status.at[i, project] == '-':
                            status.at[i, project] = 'ready'

                else:
                    status.at[step, project] = f'failed {step_status}'
            exit_log.close()

            exit_log = open(PROJECTS_PATH / project / 'exit_log.txt', 'w')
            exit_log.write('')
            exit_log.close()

            status.to_csv(status_path, sep='\t')
        except:
            #  print('exit log open exception')
            ''
        data = pd.concat([data, status], axis=1)
    return data


def draw_table(data):
    """Make a string containing the status table to print.
    Get status dataframe."""

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
    """Update status and print the interface."""

    data = update_status()
    table = draw_table(data)

    print('Database:', database)
    print(table)
    print(PROMPT)
    # for i in log.history:
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
    if mode == 'default':
        script = SCRIPTS[step]
        path = SCRIPTS_PATH / script
    elif mode == 'cluster':
        script = BATCH_SCRIPTS[step]
        path = BATCH_SCRIPTS_PATH / script

    status_path = PROJECTS_PATH / project / 'status.txt'

    print(path)
    extension = str(path).split('.')[-1]
    if extension == 'batch':
        program = 'sbatch'
    elif extension == 'py':
        program = 'python3'
    elif extension == 'sh':
        program = 'bash'

    process = subprocess.run([program, path, project, database] + list(args))  # , capture_output=True

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
        if len(args) < n_args:  # TODO: it should check if function fits a preset argument range
            log.out(f'{n_args} arguments expected, {len(args)} passed')
        else:
            action(*args)

    return wrapper


@call
def start_all_steps():
    """Start steps that have status `ready` in all projects."""
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
    """Start `step` of a `project` with specified additional arguments."""

    log.out(f'Starting Step {step} on {project}')
    start_process(int(step), project, database, *args)
    log.out('Started')


@call
def create_project(project, inp_path):
    """Create new project named `project` using query specified in `inp_path`"""

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
    """Delete specified project."""

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


# Whether the mode is default or "cluster"
# In default mode the scripts are run directly
# In cluster mode batch scripts are run
mode = 'default'

# constants
# paths
PROJECTS_PATH = Path("projects/")
BATCH_SCRIPTS_PATH = Path("batch_scripts/")
SCRIPTS_PATH = Path("scripts/")

# Scripts to be run in the default mode
SCRIPTS = {
    1: '01_blastp.sh',
    2: '02_mk_blastp_df.py',
    3: '03_plot_blastp.py',
    4: '04_filter_blastp.py',
    5: '05_cdhit.sh',
    6: '06_mk_gen_consistent.py',
    7: '07_cluster_dict.py',
    8: '08_align.sh',
    9: '09_get_translations.py',
    10: '10_get_gen_context.py',
    11: '11_trim.sh',
    12: '12_hmmscan.sh',
    13: '13_prottest.sh',
    14: '14_get_domains.py',
    15: '15_raxml.sh'
}

# Scripts to be run in the cluster mode
BATCH_SCRIPTS = {
    1: '01_blastp.batch',
    2: '02_mk_blastp_df.batch',
    3: '03_plot_blastp.batch',
    4: '04_filter_blastp.batch',
    5: '05_cdhit.batch',
    6: '06_mk_gen_consistent.batch',
    7: '07_cluster_dict.batch',
    8: '08_align.batch',
    9: '09_get_translations.batch',
    10: '10_get_gen_context.batch',
    11: '11_trim.batch',
    12: '12_hmmscan.batch',
    13: '13_prottest.batch',
    14: '14_get_domains.batch',
    15: '15_raxml.batch'
}

# Prompt listing available commands
PROMPT = '''
Prompt:
"n [project_name] [query_path]": create new project [project_name] with query at [query_path]
"s [project] [step]": start [project] [step]
"a": start all ready steps in all projects
"d [project]": delete [project]
"q": quit
'''

# This is written into the project status file `status.txt` when project is created
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

# Steps of the analysis. Used to create status dataframe
STEPS = {'Step': ['1. Blast',
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

# Dependencies of the steps on each other
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
database = 'family'  # this is hardcoded now! should be specified by user
folders = [name for name in os.listdir(PROJECTS_PATH) if os.path.isdir(os.path.join(PROJECTS_PATH, name))]
projects = [i for i in folders if i[0] != '.']

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
