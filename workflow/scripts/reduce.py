""""""

import shutil
from pathlib import Path

# Snakemake parameters
PWMS = snakemake.input[0]  # type: ignore
TFBS = snakemake.input[1]  # type: ignore
PROFILES = snakemake.input[2]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def main() -> None:
    """Main program"""
    with open(PROFILES, "r") as f:
        for line in f:
            tf_name, profile = line.split()[:2]
            for data_type in ["pwms", "tfbs"]:
                # Define the source and destination directories
                src_dir = Path(f"{snakemake.input[data_type]}/{tf_name}/{profile}")  # type: ignore
                dst_dir = Path(f"{snakemake.output[data_type]}/{tf_name}/{profile}")  # type: ignore
                # Create the destination directory if it doesn't exist
                dst_dir.mkdir(parents=True, exist_ok=True)
                # Copy each file in the source directory to the destination directory
                for src_file in src_dir.iterdir():
                    if src_file.is_file():
                        shutil.copy2(src_file, dst_dir / src_file.name)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
