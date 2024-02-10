""""""

import shutil
import tarfile
from pathlib import Path
import pandas as pd


# Snakemake parameters
IP_TAR = snakemake.input[0] # type: ignore
IP_MAP = snakemake.input[1] # type: ignore
OP_DIR = snakemake.output[0] # type: ignore


# ------------- #
# Functions     #
# ------------- #

def read_biosample_ids(filepath: str) -> pd.DataFrame:
    """Reads mapping file and returns dataframe."""
    return pd.read_csv(filepath, sep="\t", usecols=["tf_id"])

def flatten_directory(ip_dir: str) -> None:
    """Flatten a directory structure to a single directory"""
    # Set up Path on input dir
    ip_folder = Path(ip_dir)

    # Iterate over all files in the directory and its subdirectories
    for src_file in ip_folder.rglob('*'):
        if src_file.is_file():
            
            # Define the destination file path
            dst_file = ip_folder / src_file.name
            
            # Move the file if it's not already in the root directory
            if src_file != dst_file:
                shutil.move(str(src_file), str(dst_file))
                
    # Now remove the empty directories
    for folder in ip_folder.rglob('*'):
        if folder.is_dir():
            folder.rmdir()

def main() -> None:
    """Main program"""
    # Biosamples IDs to extract
    biosamples = read_biosample_ids(IP_MAP)
 
    # Ensure OP directory exists
    Path(OP_DIR).parent.mkdir(parents=True, exist_ok=True)
    
    # Extract all files from tar archive
    with tarfile.open(IP_TAR, "r") as tar:
        # Extract only files with biosample ids
        biosample_ids = biosamples["tf_id"].tolist()
        for member in tar.getmembers():
            if any(biosample_id in member.name for biosample_id in biosample_ids) and member.isfile():
                tar.extract(member, path=Path(OP_DIR).parent)
    
    # Flatten directory
    flatten_directory(OP_DIR)

# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
