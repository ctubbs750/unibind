""""""

from pandas import read_csv

# Snakemake parameters
PROFILES_MAP = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def save_mapping(mapping, output):
    """Save the final mapping to a TSV file."""
    mapping.to_csv(output, index=False, sep="\t", header=False)


def main() -> None:
    """
    Main function to read the profiles map and extract profile version.
    """
    # Define column names and types
    names = [
        "tf_name",
        "profile",
        "species",
        "class",
        "family",
        "source",
        "profile_length",
    ]
    dtype = {name: str for name in names}

    # Read the profiles map
    mapping = read_csv(
        PROFILES_MAP,
        header=None,
        sep="\t",
        names=names,
        dtype=dtype,
    )

    # Extract profile version
    mapping["version"] = mapping["profile"].str.split(".", expand=True)[1]

    # Extract linked names [id1]::[id2]
    mapping["tf_name"] = mapping["tf_name"].str.split(":", expand=True)[0]

    # A handful of manual corrections
    corrections = {
        "TBXT": "T",
        "PPARA": "RXRA",
        "NR4A2": "RXRA",
        "NR1H4": "RXRA",
        "FOSB": "JUN",
    }
    mapping["tf_name"].replace(corrections, inplace=True)

    # Reduce set to human
    mapping = mapping[mapping["species"] == "Homo sapiens"]

    # Add a leading 1 to source for sorting
    mapping.loc[mapping["source"] == "ReMap", "source"] = "1ReMap"
    mapping.loc[mapping["source"] == "ENCODE", "source"] = "1ENCODE"

    # Sort by most recent version first, then source of ReMap over individual experiment
    mapping.sort_values(
        by=["tf_name", "version", "source"], ascending=[True, False, True], inplace=True
    )

    # Drop duplicates to get final set
    mapping.drop_duplicates(subset=["tf_name"], keep="first", inplace=True)

    # Save final mapping
    save_mapping(mapping, OUTPUT)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
