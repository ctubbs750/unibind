""""""

from pathlib import Path
import os
import pandas as pd

# Snakemake parameters
IP_FOLDER = snakemake.input[0]  # type: ignore
OP_FOLDER = snakemake.output[0]  # type: ignore
EXTENSION = snakemake.params[0]  # type: ignore

# ------------- #
# Functions     #
# ------------- #


def organize_files(extension: str, ip_folder: str, op_folder: str) -> None:
    """Organizes Unibind download directories into TF/profile/datasets format"""
    # Get all files with the given extension
    paths = list(Path(ip_folder).rglob(f"*.{extension}"))
    
    print("Dflsak")
    print(ip_folder, op_folder, paths)

    # Create a DataFrame with file paths and info
    df = pd.DataFrame({"path": paths, "info": [p.stem for p in paths]})

    # Split the info into different fields
    df[["exp_id", "celltype", "tf_name", "profile_root", "profile_stem", "damo"]] = df[
        "info"].str.split(".", expand=True)

    # Create a profile tag
    df["profile"] = df["profile_root"] + "." + df["profile_stem"]

    # Group by all unique TF/profiles
    grouped = df.groupby(["tf_name", "profile"])["path"].apply(list).reset_index()

    # Loop over groups and move files
    for row in grouped.itertuples():
        dest_dir = Path(op_folder) / row.tf_name / row.profile
        dest_dir.mkdir(parents=True, exist_ok=True)

        for path in row.path:
            os.rename(path, dest_dir / path.name)


def main():
    """Main program"""
    organize_files(EXTENSION, IP_FOLDER, OP_FOLDER)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
