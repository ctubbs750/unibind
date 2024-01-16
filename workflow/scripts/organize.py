#!/usr/bin/python

from pathlib import Path
from shutil import move
from pandas import DataFrame

# Snakemake parameters
IP_FOLDER = snakemake.input[0]  # type: ignore
OP_FOLDER = snakemake.output[0]  # type: ignore

###
# Functions
###


def main():
    """"""
    # File imfo
    paths = [i for i in Path(IP_FOLDER).glob("*/*.bed")]
    stems = [i.stem for i in paths]

    # Convert to frame
    df = DataFrame(list(zip(paths, stems)), columns=["path", "info"])

    # Break UniBind label into fields
    df[["exp_id", "celltype", "tf_name", "profile_root", "profile_stem", "damo"]] = df[
        "info"
    ].str.split(".", expand=True)

    # Make profile tag
    df["profile"] = df["profile_root"] + "." + df["profile_stem"]

    # Group by all unique TF/profiles
    grouped = df.groupby(["tf_name", "profile"])["path"].unique().reset_index()

    # Loop over groups
    for row in grouped.itertuples():
        tf_name = row.tf_name
        profile = row.profile
        filepaths = row.path

        # Move files
        Path(f"{OP_FOLDER}/{tf_name}/{profile}").mkdir(parents=True, exist_ok=True)
        for path in filepaths:
            filename = path.stem
            move(path, f"{OP_FOLDER}/{tf_name}/{profile}/{filename}")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
