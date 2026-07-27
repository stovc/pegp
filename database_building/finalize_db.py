"""Concatenate database fragments and remove their source folders."""

from argparse import ArgumentParser
from pathlib import Path

from utils import concat_files_from_folder, delete_folder


FOLDERS_TO_CONCATENATE_CSV = ["annotation", "sequence", "translation"]


def parse_arguments():
    parser = ArgumentParser(
        description="Concatenate database fragments and remove their source folders."
    )
    parser.add_argument(
        "database",
        help="Path to the database directory containing the fragment folders.",
    )
    return parser.parse_args()


def finalize_database(database_path):
    """Concatenate database fragments and delete the fragment folders."""
    concat_files_from_folder(
        folder=database_path / "protein",
        extension="faa",
    )

    for folder in FOLDERS_TO_CONCATENATE_CSV:
        concat_files_from_folder(
            folder=database_path / folder,
            extension="csv",
        )

    for folder in FOLDERS_TO_CONCATENATE_CSV + ["protein"]:
        delete_folder(database_path / folder)


if __name__ == "__main__":
    args = parse_arguments()
    finalize_database(Path(args.database))
