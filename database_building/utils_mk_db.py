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


def circular_slice(s, left, right):
    """Slice a string as circular."""
    if left < 0:
        return s[left:] + s[:right]
    elif right > len(s):
        return s[left:] + s[:right - len(s)]
    else:
        return s[left:right]


def get_first(dict_arg, key):
    """Return first element of a dict item if possible. If not subscriptable, return the item."""
    get = dict_arg.get(key)
    try:
        return get[0]
    except:
        return get


def concatenate(folder, extension):
    """Concatenates all files of 'folder' into a single file of 'extension' extension."""
    print('folder', folder)
    record_list = os.listdir(database_path / folder)
    record_list.sort()

    with open(database_path / f'{folder}.{extension}', 'w') as outfile:
        for f in record_list:
            with open(database_path / folder / f, 'r') as infile:
                outfile.write(infile.read())


def concatenate_csv(folder):
    """Concatenates all files of 'folder' into a single csv file."""
    print('folder', folder)
    record_list = os.listdir(database_path / folder)
    record_list.sort()

    df_concat = pd.concat([pd.read_csv(database_path / folder / f) for f in record_list], ignore_index=True)
    df_concat.to_csv(database_path / f'{folder}.csv', index=False)