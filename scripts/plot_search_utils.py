from Bio import SeqIO
import constants as c


def generate_colormap(list_of_values):
    """Generate a colormap dictionary mapping each value to a corresponding color."""
    return {value: c.SCATTERPLOT_COLORS[i] for i, value in enumerate(list_of_values)}


def put_other_to_the_end(data, distribution):
    """Put 'Other' to the end of `data` if present
    data - array with features of hits
    distribution - number of hits with particular feature
    """

    if 'Other' in data:
        other_position = data.index('Other')
        other_count = distribution[other_position]
        del data[other_position]
        del distribution[other_position]
        data.append('Other')
        distribution.append(other_count)


def put_other_to_end(value_counts):

    index_list = value_counts.index.tolist()
    sorted_index_list = sorted(index_list)

    # Check if 'Other' is in the list
    if 'Other' in sorted_index_list:
        # Remove 'Other' from the list
        sorted_index_list.remove('Other')
        # Append 'Other' to the end of the list
        sorted_index_list.append('Other')

    sorted_value_counts = value_counts[sorted_index_list]

    return sorted_value_counts


def sort_list_with_other_at_end(input_list):
    sorted_list = sorted(input_list)

    # Check if 'Other' is in the list
    if 'Other' in sorted_list:
        # Remove 'Other' from the list
        sorted_list.remove('Other')
        # Append 'Other' to the end of the list
        sorted_list.append('Other')

    return sorted_list


def count_fasta_records(file_path):
    return sum(1 for _ in SeqIO.parse(file_path, "fasta"))


