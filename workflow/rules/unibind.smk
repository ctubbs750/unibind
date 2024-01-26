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


rule download_unibind_pwms:
    message:
        "Downloads all PWMS from UniBind database"
    output:
        PWMS_DOWNLOAD,
    params:
        url=PWMS_URL,
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="download_unibind_pwms"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="download_unibind_pwms"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        "curl -o {output} {params.url}"


rule unpack_unibind_pwms:
    message:
        "Unpacks UniBind PWM download"
    input:
        rules.download_unibind_pwms.output,
    output:
        temp(directory(PWMS_UNPACKED)),
    params:
        outdir=lambda w, input: os.path.dirname(input[0]),
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="unpack_unibind_pwms"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="unpack_unibind_pwms"),
    conda:
        "../envs/unibind.yaml"
    shell:
        "tar -xzf {input} -C {params.outdir}"


rule organize_pwms:
    message:
        "Organizes UniBind damo PWMS by TF and profile"
    input:
        rules.unpack_unibind_pwms.output,
    output:
        temp(directory(PWMS_ORGANIZED)),
    params:
        extension=EXTENSIONS["pwms"],
        script="../scripts/organize.py",
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="organize_pwms"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="organize_pwms"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "{params.script}"


rule download_unibind_tfbs:
    message:
        "Downloads all TFBSs from UniBind database"
    output:
        TFBS_DOWNLOAD,
    params:
        url=TFBS_URL,
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="download_unibind_tfbs"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="download_unibind_tfbs"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        "curl -o {output} {params.url}"


rule unpack_unibind_tfbs:
    message:
        "Unpacks UniBind TFBS download"
    input:
        rules.download_unibind_tfbs.output,
    output:
        temp(directory(TFBS_UNPACKED)),
    params:
        outdir=lambda w, input: os.path.dirname(input[0]),
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="unpack_unibind_tfbs"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="unpack_unibind_tfbs"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        "tar -xzf {input} -C {params.outdir}"


rule organize_tfbs:
    message:
        "Organizes UniBind damos by TF and profile"
    input:
        rules.unpack_unibind_tfbs.output,
    output:
        temp(directory(TFBS_ORGANIZED)),
    params:
        extension=EXTENSIONS["tfbs"],
        script="../scripts/organize.py",
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="organize_tfbs"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="organize_tfbs"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "{params.script}"


rule profile_mapping:
    message:
        "Confirms each profile is homo sapiens."
    input:
        rules.organize_pwms.output,
    output:
        PROFILE_MAPPING,
    params:
        script="../scripts/mapping.py",
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="profile_mapping"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="profile_mapping"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "{params.script}"


rule filter_tf_targets:
    message:
        "Reduces list of profiles to just those that are human and non-redundant."
    input:
        rules.profile_mapping.output,
    output:
        PROFILE_MAPPING_FILTERED,
    params:
        script="../scripts/profiles.py",
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="filter_tf_targets"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="filter_tf_targets"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "{params.script}"


rule reduce_data_to_filter:
    message:
        "Reduces downloaded data to those that are human"
    input:
        pwms=rules.organize_pwms.output,
        tfbs=rules.organize_tfbs.output,
        profiles=rules.filter_tf_targets.output,
    output:
        pwms=directory(PWMS_ORGANIZED_FILTER),
        tfbs=directory(TFBS_ORGANIZED_FILTER),
    params:
        script="../scripts/reduce.py",
    log:
        stdout=expand("workflow/logs/{rule}.stdout", rule="reduce_data_to_filter"),
        stderr=expand("workflow/logs/{rule}.stderr", rule="reduce_data_to_filter"),
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "{params.script}"
