import logging

import pandas as pd

import itertools
import os
import platform
import subprocess


def distance_between_sequences(start1, end1, start2, end2, circular_length=None):
    """
    Return distance between two SeqFeatures. Supports circular nucleic acids.
    Nucleic acid is supposed to be linear unless circular_length argument is passed.

    Parameters
    ----------
    start1, end1 -- start and end nucleotides of the 1st SeqFeature
    start2, end2 -- start and end nucleotides of the 2nd SeqFeature
    circular_length -- length of a circular nucleic acid; nucleic acid supposed linear if None

    Returns
    -------
    distance -- distance between two SeqFeatures.
    """
    if circular_length is None:
        x = min(abs(start1 - end2),
                abs(end1 - start2))
    else:
        x = min(abs(start1 - end2),
                abs(end1 - start2),
                abs(circular_length - end1 + start2),
                abs(circular_length - end2 + start1))
    return x


def circular_slice(string, left, right):
    """Slice a string as circular."""
    if left < 0:
        return string[left:] + string[:right]
    elif right > len(string):
        return string[left:] + string[:right - len(string)]
    else:
        return string[left:right]


def get_first(dict_arg, key):
    """Return first element of a dict item if possible. If not subscriptable, return the item."""
    get = dict_arg.get(key)
    try:
        return get[0]
    except:
        return get


def concat_files_from_folder(folder, extension):
    """Concatenates all files of 'folder' into a single file of 'extension' extension."""
    logging.info('concatenating folder', folder)
    record_list = os.listdir(folder)
    record_list.sort()

    if extension == 'csv':
        df_concat = pd.concat([pd.read_csv(folder / f, engine='python') for f in record_list], ignore_index=True)
        df_concat.to_csv(f'{folder}.csv', index=False)

    else:
        with open(f'{folder}.{extension}', 'w') as outfile:
            for f in record_list:
                with open(folder / f, 'r') as infile:
                    outfile.write(infile.read())

    return None


def delete_folder(folder_path):
    try:
        if platform.system() == "Windows":
            subprocess.run(["rmdir", "/s", "/q", folder_path], shell=True, check=True)
        else:
            subprocess.run(["rm", "-rf", folder_path], check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"Error: {e}")
    except Exception as e:
        logging.error(f"An unexpected error occurred: {e}")


def make_id_iterator(prefix, symbols, suffix_length):

    combinations = itertools.product(symbols, repeat=suffix_length)

    for combination in combinations:
        yield prefix + ''.join(combination)
