""""""

from pandas import DataFrame, read_csv
from numpy import int64
from pathlib import Path
from scipy.stats import zscore


# Snakemake parameters
PWMS_DIR = snakemake.input[0]  # type: ignore
TFBS_DIR = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


# ------------- #
# Functions     #
# ------------- #

def pwm_to_IntLogOdds(filepath: str) -> DataFrame:
    """Converts PWM to integer log-odds format"""
    # Read PWM, values are assumed to be log-odds
    pwm = read_csv(filepath, header=None, sep="\t")
    # Round to integer, transpose for PWMScan format
    pwm = pwm.round().astype(int64).T
    return pwm
    

def calculate_ic(filepath: str) -> tuple:
    """Approximates IC from masked PWM"""
    pwm_to_convert = pwm_to_IntLogOdds(filepath)
    # length of pwm
    pwm_len = pwm_to_convert.shape[0]
    #
    num_pos = len([i for i in pwm_to_convert.max(axis=1) if i >0])
    # proportion of positions with positive values
    prop = num_pos / pwm_len
    
    return pwm_to_convert.max(axis=1).mean(), prop

def main():
    # Get all biosamples linked to profile
    biosamples = [i.name for i in Path(PWMS_DIR).glob("*/*/*")]

    # Creat ampping of profiles to biosamples
    profile_map = {}
    for biosample in biosamples:
        tf_name = biosample.split(".")[-5]
        profile = ".".join(biosample.split(".")[-4:-2])
        # Get IC for each profile
        ic, prop_pos = calculate_ic(f"{PWMS_DIR}/{tf_name}/{profile}/{biosample}")
        #
        if profile not in profile_map:
            profile_map[profile] = []
            profile_map[profile].append( [biosample, ic, prop_pos])
        else:
            profile_map[profile].append( [biosample, ic, prop_pos])
    
    #
    fails = []
    for profile in profile_map:
        # Get biosamples
        biosamples = [i[0] for i in profile_map[profile]]
        # Tf name
        tf_name = biosamples[0].split(".")[-5]
        # Get all ICs for profile
        ics = [float(i[1]) for i in profile_map[profile]]
        prop_pos = [float(i[2]) for i in profile_map[profile]]
        # Convert to zscores
        z_ics = zscore(ics)
        for index, entry in enumerate(zip(z_ics, prop_pos)):
            z = entry[0]
            prop_pos = entry[1]
            if (z < -2 or prop_pos<1):
                # Now clean out profile and TFBS dirs
                pwms_path = Path(f"{PWMS_DIR}/{tf_name}/{profile}/{biosamples[index]}")
                tfbs_path = Path(f"{TFBS_DIR}/{tf_name}/{profile}/{biosamples[index][:-4]}.bed")
                # Remove files
                pwms_path.unlink(missing_ok=True)
                tfbs_path.unlink(missing_ok=True)
                # Save fail
                fails.append([profile, biosamples[index], z, prop_pos])
                
    # Make fails into matrix
    fails = DataFrame(fails)
    fails.columns = ["profile", "biosample", "z_ic", "prop_pos"]

    # Save fails
    fails.to_csv(OUTPUT, index=False, sep="\t")
    
    # As a last step, remove empty directories from TFBS dirs and PWMS dirs
    for tf in Path(TFBS_DIR).glob("*"):
        for profile in tf.glob("*"):
            if len(list(profile.glob("*"))) == 0:
                profile.rmdir()
                tf.rmdir()
        

# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()