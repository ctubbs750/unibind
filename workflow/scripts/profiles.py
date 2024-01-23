#!/usr/bin/python

from pandas import read_csv

# Snakemake parameters
PROFILES_MAP = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore

###
# Functions
###


def main() -> None:
    """"""
    # Read map
    names = ["tf_name", "profile", "species", "class", "family", "source"]
    dtype = [str, str, str, str, str, str]
    mapping = read_csv(
        PROFILES_MAP,
        header=None,
        delim_whitespace=True,
        names=names,
        dtype=dict(zip(names, dtype)),
    )

    # Extract profile version
    mapping["version"] = [i[1] for i in mapping["profile"].str.split(".")]

    #  Extract linked names [id1]::[id2]
    mapping["tf_name"] = [i[0] for i in mapping["tf_name"].astype(str).str.split(":")]

    # A handful of manual corrections
    # TBXT=T, PPARA=RXRA keep RXRA, NR4A2==RXRA, NR1H4=RXRA, FOSB==JUN
    mapping["tf_name"] = mapping["tf_name"].replace(
        {"TBXT": "T", "PPARA": "RXRA", "NR4A2": "RXRA", "NR1H4": "RXRA", "FOSB": "JUN"}
    )

    # Reduce set to human
    mapping = mapping[mapping["species"] == "Homo sapiens"]

    # Add a leading 1 to source for sorting
    mapping["source"] = ["1" + i if i == "ReMap" else i for i in mapping["source"]]

    # Selects on most recent version first, then source of ReMap over individual experiment
    mapping.sort_values(
        by=["tf_name", "version", "source"], ascending=[True, False, True], inplace=True
    )

    # Drop duplicates to get final set
    mapping.drop_duplicates(subset=["tf_name"], keep="first", inplace=True)

    # Save final mapping
    mapping.to_csv(OUTPUT, index=False, sep=" ", header=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
