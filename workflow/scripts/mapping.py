""""""

import tarfile
import coreapi
from pathlib import Path
import pandas as pd


# Snakemake parameters
IP_TAR = snakemake.input[0] # type: ignore
OUTPUT = snakemake.output[0] # type: ignore


# ------------- #
# Functions     #
# ------------- #


def get_biosamples(tar_path: str) -> None:
    """Returns list of biosamples in tar archives."""
    biosamples = []
    # Open tar archive
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            # Check if member is a file and ends with pwm
            if member.isfile() and member.name.endswith("pwm"):
                biosamples.append(Path(str(member)).parts[1])
    # Return list of biosamples
    return biosamples


def extract_unibind_info(data: dict) -> dict:
    """Extracts info from api pulldown"""
    # Parse results from data query
    return  {
        "tf_id": data.get("tf_id", "NaN"),
        "tf_name":  data.get("tf_name", "NaN"),
        "biological_condition": data.get("biological_condition", "NaN"),
        "jaspar_id": data.get("jaspar_id", "NaN")[0],
        "jaspar_version": int(data.get("tfbs")[0]["DAMO"][0]["jaspar_version"]), #type: ignore
        "score_threshold": float(data.get("tfbs")[0]["DAMO"][0]["score_threshold"]) ,#type: ignore
        "total_tfbs": int(data.get("tfbs")[0]["DAMO"][0]["total_tfbs"]),#type: ignore
        "total_peaks": data.get("total_peaks", "NaN")
    }
    
def extract_jaspar_info(data: dict) -> dict:
    """Extracts jaspar info from api pulldown"""
    return {
        "name": data.get("name", "NaN"),
        "species": data.get("species", [{"name": "NaN"}])[0]["name"],
        "class": data.get("class", ["NaN"])[0],
        "family": data.get("family", ["NaN"])[0],
        "source": data.get("source", "NaN"),
        "length": len(data.get("pfm", {"A": []})["A"]),
    }


def main() -> None:
    """Main program"""
    # Get all file names from input dir
    biosamples = get_biosamples(IP_TAR)

    # Initialize a client & load the schema document
    client = coreapi.Client()
    schema = client.get("https://unibind.uio.no/api/v1/docs")
    action = ["datasets", "read"]

    # Loop over biosamples to get info
    biosample_info = []
    for biosample in biosamples:
        # Extract the biosample ID and cell type
        biosample_id = ".".join(biosample.split(".")[:3])
        biosample_ct = biosample.split(".")[1]
        print(biosample_id)
        # Interact with the API endpoint
        params = {"tf_id": biosample_id}
        result = client.action(schema, action, params=params)
        # Extract the data
        info = extract_unibind_info(result)
        info["biosample_ct"] = biosample_ct
        biosample_info.append(info)
        
    # Convert to dataframe
    bi = pd.DataFrame(biosample_info)

    # Convert empty lists in biological condiations to NaN
    bi["biological_condition"] = bi["biological_condition"].apply(lambda x: "NaN" if len(x) == 0 else x)

    # Get list of all unique profile IDs
    profile_ids = list(bi["jaspar_id"].unique())

    # Initialize JASPAR client & load the schema document
    client = coreapi.Client()
    schema = client.get("https://jaspar.elixir.no/api/v1/docs")
    action = ["matrix", "read"]

    # Interact with the API endpoint
    profiles_info = []
    for profile in profile_ids:
        print(profile)
        # Pull profile info
        params = {"matrix_id": profile}
        # Extract the data
        try:
            result = client.action(schema, action, params=params)
            info = extract_jaspar_info(result)
        except:
            info = {"name": "NaN",
                    "species": "NaN",
                    "class": "NaN",
                    "family": "NaN",
                    "source": "NaN",
                    "length": "NaN"}
        info["profile"] = profile
        profiles_info.append(info)
        
    # Convert to dataframe
    pi = pd.DataFrame(profiles_info)

    # Merge bi and pi on profile
    mapping = pd.merge(bi, pi, left_on="jaspar_id", right_on="profile", how="left")
    
    # Reset index
    mapping.reset_index(drop=True, inplace=True)
    
    # Save mapping
    mapping.to_csv(OUTPUT, index=False, sep="\t")


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
