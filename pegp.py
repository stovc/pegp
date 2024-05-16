import os
import subprocess
import pandas as pd
from prettytable import PrettyTable
from pathlib import Path
from colorama import Fore, Style

from config import *


# INPUT OUTPUT
def update_status():
    """Update the status of each step of each project.
    Iterate project directories, look at `status.txt` files.
    Return dataframe with status of steps in projects."""

    status_data = pd.DataFrame(STEPS, index=list(range(1, 18)))

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
            print('exit log open exception')

        status_data = pd.concat([status_data, status], axis=1)
    return status_data


def draw_status_table(status_data):
    """Print a string containing the status table to print.
    Get status dataframe."""

    table = PrettyTable(['Step'] + projects)
    table.align["Step"] = "l"
    for row in status_data.iterrows():
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
    print(table)


def update_screen():
    """Update status and print the interface."""
    print('Database:', database)
    data = update_status()
    draw_status_table(data)
    print(PROMPT)
    # for i in log.history:
    #    print(i)


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
            log.save_and_print(f'{n_args} arguments expected, {len(args)} passed')
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
        log.save_and_print(f'Starting Step {step} on {project}')
        start_process(step, project, database)
        log.save_and_print('Started')


@call
def start_step(project, step, *args):
    """Start `step` of a `project` with specified additional arguments."""

    log.save_and_print(f'Starting Step {step} on {project}')
    start_process(int(step), project, database, *args)
    log.save_and_print('Started')


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

        log.save_and_print(f'Created a new project {project}')
    else:
        log.save_and_print(f'Project with name {project} already exists')


@call
def delete_project(project):
    """Delete specified project."""

    delete_project_path = PROJECTS_PATH / project
    if os.path.exists(delete_project_path):
        os.system(f'rm -R {delete_project_path}')

        projects.remove(project)

        update_screen()

        log.save_and_print(f'Deleted the project {project}')
    else:
        log.save_and_print(f'Project with name {project} does not exist')


@call
def update():
    update_screen()


@call
def quitx():
    quit()


class Log:
    """Log containing the history of the input and output of the session which is printed with screen update."""
    history = []

    def save(self, message):
        self.history.append(f'>>> {message}')

    def save_and_print(self, message):
        print(message)
        self.history.append(message)


# Configs

# Whether the mode is default or "cluster"
# In default mode the scripts are run directly
# In cluster mode batch scripts are run
mode = DEFAULT_MODE

# database for homology search with metadata
database = DEFAULT_DATABASE  # this is hardcoded now! should be specified by user


# Init

folders = [name for name in os.listdir(PROJECTS_PATH) if os.path.isdir(os.path.join(PROJECTS_PATH, name))]
projects = [i for i in folders if i[0] != '.']

command = ''
log = Log()

update_screen()

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
    log.save(command)
    command_args = command.split(' ')
    if command_args[0] in commands.keys():
        action = commands[command_args[0]]
        action(command)
    else:
        log.save_and_print(f'unknown command "{command}"')
