""""""

import pandas as pd
from scipy.stats import zscore


# Snakemake parameters
IP_MAP = snakemake.input[0] # type: ignore
OUTPUT = snakemake.output[0] # type: ignore
REPORT = snakemake.output[1] # type: ignore


# ------------- #
# Functions     #
# ------------- #


def read_mapping(filepath: str) -> pd.DataFrame:
    """Reads mapping file and returns dataframe."""
    return pd.read_csv(filepath, sep="\t")


def main() -> None:
    """Main program"""
    # read in mapping
    mapping = read_mapping(IP_MAP)
    
    # NOTE: Not sure where these duplicate rows are being introduced
    mapping = mapping.drop_duplicates(keep="first")      
    
    # Get dict of all most recent profile versions per TF
    profile_map = mapping.groupby("tf_name")["profile"].max().to_dict()

    # Z-score for total_peaks within tf_name group
    mapping['total_peaks_zscore'] = mapping.groupby("tf_name", group_keys=False)["total_peaks"].apply(zscore)

    # Identify all biosamples with a bioloical condition (treatment)
    mapping['biological_condition_pass'] = [1 if x == "nan" else 0 for x in mapping['biological_condition'].astype(str)]

    # Filter to profiles derived from species homo sapiends
    mapping["species_pass"] = [1 if x=="Homo sapiens" else 0 for x in mapping["species"]]

    # Filter to most recent jaspar profile
    mapping["profile_pass"] = [1 if profile==profile_map[tf_name] else 0 for tf_name, profile in zip(mapping["tf_name"], mapping["profile"])]

    # Filter to tfbs peak set abs z score < 2 AND have at least 1000 peaks
    mapping['tfbs_peaks_pass'] = [1 if abs(zscore) < 2 and count >=1000 else 0 for zscore, count in zip(mapping["total_peaks_zscore"], mapping["total_peaks"])]

    # Now reduce to where all filters are pass
    filters = ["biological_condition_pass", "species_pass", "profile_pass", "tfbs_peaks_pass"]
    filtered = mapping.loc[(mapping[filters[0]]==1) & (mapping[filters[1]]==1) & (mapping[filters[2]]==1)& (mapping[filters[3]]==1)]
    
    # Save filtered
    filtered.to_csv(OUTPUT, index=False, sep="\t")

    # Make report
    with open(REPORT, 'w') as f:
        # Total number of samples
        print(f"Total number of biosamples: {len(mapping)}", file=f)
        for criteria in filters:
            print(f"{criteria} FAILS: {len(mapping[mapping[criteria]==0])}", file=f)
        # total remainging
        print(f"Total number of biosamples remaining: {len(filtered)}", file=f)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
