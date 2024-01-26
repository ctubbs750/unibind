""""""

import re
from pathlib import Path
import coreapi

# Snakemake parameters
IP_DIR = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #


def get_profile_info(profile: str) -> dict:
    """"""
    # Initialize a client & load the schema document
    client = coreapi.Client()
    schema = client.get("https://jaspar.elixir.no/api/v1/docs")

    # Interact with the API endpoint
    action = ["matrix", "read"]
    params = {"matrix_id": profile}

    # Try to get the data
    try:
        data = client.action(schema, action, params=params)
    except:
        data = {}

    # Return the desired info
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
    # Define the root directory and the desired depth
    root = Path(IP_DIR)
    regex = re.compile(r"MA\d+\.\d+")

    # Get list off all profiles
    profiles = [
        str(d.name) for d in root.rglob("*") if d.is_dir() and regex.match(d.name)
    ]

    # Get profile info
    with open(OUTPUT, "w") as f:
        for profile in profiles:
            # Get profile info
            info = get_profile_info(profile)
            f.write(
                f"{info['name']}\t{profile}\t{info['species']}\t{info['class']}\t{info['family']}\t{info['source']}\t{info['length']}\n"
            )


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
