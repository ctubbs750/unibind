from os import listdir, path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #


PWMS_URL = config["urls"]["pwms"]
TFBS_URL = config["urls"]["tfbs"]
INSTALL_DIR = config["install_dir"]
PROCESS_DIR = config["process_dir"]

# ------------- #
# I/O           #
# ------------- #

# Raw PWM and TFBS download
PWMS_DOWNLOAD = os.path.join(INSTALL_DIR, "damo_hg38_PWMs.tar.gz")
TFBS_DOWNLOAD = os.path.join(INSTALL_DIR, "damo_hg38_TFBS_per_TF.tar.gz")

# Unpacked PWMs
PWMS_UNPACKED = os.path.join(INSTALL_DIR, "damo_hg38_PWMs")
TFBS_UNPACKED = os.path.join(INSTALL_DIR, "damo_hg38_TFBS_per_TF")

# Orgnaized PWMs
PWMS_ORGANIZED = os.path.join(PROCESS_DIR, "damo_hg38_PWMS")
TFBS_ORGANIZED = os.path.join(PROCESS_DIR, "damo_hg38_TFBS")

# All species IDs for each profile and just those that are human
PROFILE_MAPPING = os.path.join(PROCESS_DIR, "profile_mapping.txt")
PROFILE_MAPPING_FILTERED = os.path.join(PROCESS_DIR, "profile_mapping_filtered.txt")

# Final human specific datasetes
PWMS_ORGANIZED_FILTER = os.path.join(PROCESS_DIR, "damo_hg38_PWMS_filter")
TFBS_ORGANIZED_FILTER = os.path.join(PROCESS_DIR, "damo_hg38_TFBS_filter")

# ------------- #
# Params        #
# ------------- #

# Extensions on PWMs and TFBS downloads
EXTENSIONS = {"pwms": "pwm", "tfbs": "bed"}

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        PWMS_ORGANIZED,
        TFBS_ORGANIZED,
        PWMS_ORGANIZED_FILTER,
        TFBS_ORGANIZED_FILTER,
    default_target: True


rule download_unibind_pwms:
    message:
        """
        Downloads all PWMS from UniBind database
        """
    output:
        PWMS_DOWNLOAD,
    params:
        url=PWMS_URL,
    log:
        stdout="workflow/logs/download_unibind_pwms.stdout",
        stderr="workflow/logs/download_unibind_pwms.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        """
        curl -o {output} {params.url}
        """


rule unpack_unibind_pwms:
    message:
        """
        Unpacks UniBind PWM download
        """
    input:
        rules.download_unibind_pwms.output,
    output:
        temp(directory(PWMS_UNPACKED)),
    params:
        outdir="resources/data/unibind",
    log:
        stdout="workflow/logs/unpack_unibind_pwms.stdout",
        stderr="workflow/logs/unpack_unibind_pwms.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        """
        mkdir -p {output} && tar -xzf {input} -C {params.outdir}
        """


rule organize_pwms:
    message:
        """
        Orgnzies UniBind damo PWMS by TF and profile
        """
    input:
        rules.unpack_unibind_pwms.output,
    output:
        temp(directory(PWMS_ORGANIZED)),
    params:
        extension=EXTENSIONS["pwms"],
    log:
        stdout="workflow/logs/organize_pwms.stdout",
        stderr="workflow/logs/organize_pwms.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "../scripts/organize.py"


rule download_unibind_tfbs:
    message:
        """
        Downloads all TFBSs from UniBind database
        """
    output:
        TFBS_DOWNLOAD,
    params:
        url=TFBS_URL,
    log:
        stdout="workflow/logs/download_unibind_tfbs.stdout",
        stderr="workflow/logs/download_unibind_tfbs.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        """
        curl -o {output} {params.url}
        """


rule unpack_unibind_tfbs:
    message:
        """
        Unpacks UniBind TFBS download
        """
    input:
        rules.download_unibind_tfbs.output,
    output:
        temp(directory(TFBS_UNPACKED)),
    params:
        outdir=INSTALL_DIR,
    log:
        stdout="workflow/logs/unpack_unibind_tfbs.stdout",
        stderr="workflow/logs/unpack_unibind_tfbs.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        """
        mkdir -p {output} && tar -xzf {input} -C {params.outdir}
        """


rule organize_tfbs:
    message:
        """
        Orgnzies UniBind damos by TF and profile
        """
    input:
        rules.unpack_unibind_tfbs.output,
    output:
        temp(directory(TFBS_ORGANIZED)),
    params:
        extension=EXTENSIONS["tfbs"],
    log:
        stdout="workflow/logs/organize_files.stdout",
        stderr="workflow/logs/organize_files.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "../scripts/organize.py"


rule profile_mapping:
    message:
        """
        Confirms each profile is homo sapiens.
        """
    input:
        rules.organize_pwms.output,
    output:
        PROFILE_MAPPING,
    log:
        stdout="workflow/logs/check_human_profile.stdout",
        stderr="workflow/logs/check_human_profile.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        """
        coreapi get https://jaspar.elixir.no/api/v1/docs/
        readarray -t profiles < <(find {input} -mindepth 2 -maxdepth 2 -name '*M*' -type d -exec basename {{}}  \;)
        for profile in "${{profiles[@]}}";
        do
            echo "$profile"
            query=$(coreapi action matrix read -p matrix_id="$profile") || query="NaN"

            tf_species=$(echo $query | jq '.species.[].name' | tr -d "\n") || species="NaN"
            tf_name=$(echo $query | jq '.name' | tr -d "\n") || tf_name="NaN"
            tf_class=$(echo $query | jq '.class[]' | tr -d "\n") || tf_class="NaN"
            tf_family=$(echo $query | jq '.family[]' | tr -d "\n") || tf_family="NaN"
            tf_source=$(echo $query | jq '.source' | tr -d "\n") || tf_source="NaN"

            printf '%s %s %s %s %s %s\n' "$tf_name" "$profile" "$tf_species" "$tf_class" "$tf_family" "$tf_source" >> {output}
        done
        """


rule filter_tf_targets:
    message:
        """
        Reduces list of profiles to just thoe that are human and non-redundant.
        """
    input:
        rules.profile_mapping.output,
    output:
        PROFILE_MAPPING_FILTERED,
    conda:
        "../envs/unibind.yaml"
    log:
        stdout="workflow/logs/filter_tf_targets.stdout",
        stderr="workflow/logs/filter_tf_targets.stderr",
    threads: 1
    script:
        "../scripts/profiles.py"


rule reduce_data_to_filter:
    message:
        """
        Reducess downloaded data to those that are human
        """
    input:
        pwms=rules.organize_pwms.output,
        tfbs=rules.organize_tfbs.output,
        profiles=rules.filter_tf_targets.output,
    output:
        pwms=directory(PWMS_ORGANIZED_FILTER),
        tfbs=directory(TFBS_ORGANIZED_FILTER),
    conda:
        "../envs/unibind.yaml"
    log:
        stdout="workflow/logs/reduce_data_to_filter.stdout",
        stderr="workflow/logs/reduce_data_to_filter.stderr",
    threads: 1
    shell:
        """
        while read p; do
            # Extract line info
            tf_name=$(echo "$p" | awk '{{print $1}}')
            echo $tf_name
            profile=$(echo "$p" | awk '{{print $2}}')
            # Setup dirs
            mkdir -p {output.pwms}/${{tf_name}}
            mkdir -p {output.tfbs}/${{tf_name}}
            # Move relevant dir/files to output
            cp -r {input.pwms}/${{tf_name}}/${{profile}} {output.pwms}/${{tf_name}}
            cp -r {input.tfbs}/${{tf_name}}/${{profile}} {output.tfbs}/${{tf_name}}
        done < {input.profiles}
        """
