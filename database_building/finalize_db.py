"""Concatenate database fragments and remove their source folders."""

from argparse import ArgumentParser
import logging
from pathlib import Path

from utils import concat_files_from_folder, delete_folder


FOLDERS_TO_CONCATENATE_CSV = ["annotation", "sequence", "translation"]
LOGGER = logging.getLogger(__name__)


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
    folders = [("protein", "faa")]
    folders.extend((folder, "csv") for folder in FOLDERS_TO_CONCATENATE_CSV)

    total_files = 0
    for folder_name, extension in folders:
        folder = database_path / folder_name
        file_count = sum(path.is_file() for path in folder.iterdir())
        output_path = database_path / f"{folder_name}.{extension}"

        LOGGER.info(
            "Concatenating %d files from '%s' into '%s'",
            file_count,
            folder,
            output_path,
        )

        def report_progress(processed_files, total_folder_files, file_name):
            LOGGER.info(
                "[%s] Processed %d/%d files: %s",
                folder_name,
                processed_files,
                total_folder_files,
                file_name,
            )

        concat_files_from_folder(
            folder=folder,
            extension=extension,
            progress_callback=report_progress,
        )
        LOGGER.info("Finished processing '%s'", folder)
        total_files += file_count

    for folder_name, _ in folders:
        folder = database_path / folder_name
        LOGGER.info("Removing source folder '%s'", folder)
        delete_folder(folder)

    LOGGER.info(
        "Finalization complete: processed %d folders and %d files",
        len(folders),
        total_files,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_arguments()
    finalize_database(Path(args.database))
