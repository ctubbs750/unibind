from os import listdir, path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

INSTALL_DIR = config["install_dir"]
PROCESS_DIR = config["process_dir"]
UNIBIND_URLS = {
    "pwms": config["urls"]["pwms"],
    "tfbs": config["urls"]["tfbs"],
}


# ------------- #
# I/O           #
# ------------- #

# Raw PWM and TFBS download
PWMS_DOWNLOAD = path.join(INSTALL_DIR, "damo_hg38_PWMS.tar.gz")
TFBS_DOWNLOAD = path.join(INSTALL_DIR, "damo_hg38_TFBS.tar.gz")

# Biosample maps
BIOSAMPLE_MAP = os.path.join(PROCESS_DIR, "biosample_map.tsv")

# Biosample map with filters labeled
BIOSAMPLE_MAP_FILTER = os.path.join(PROCESS_DIR, "biosample_map.filter.tsv")

# Extracted data
EXTRACTED_PWMS = os.path.join(PROCESS_DIR, "damo_hg38_PWMs")
EXTRACTED_TFBS = os.path.join(PROCESS_DIR, "damo_hg38_TFBS_per_TF")

# # Unpacked PWMs
# PWMS_UNPACKED = os.path.join(INSTALL_DIR, "damo_hg38_PWMs")
# TFBS_UNPACKED = os.path.join(INSTALL_DIR, "damo_hg38_TFBS_per_TF")

# # # Orgnaized PWMs
# # PWMS_ORGANIZED = os.path.join(PROCESS_DIR, "damo_hg38_PWMS")
# # TFBS_ORGANIZED = os.path.join(PROCESS_DIR, "damo_hg38_TFBS")

# # Orgnaized PWMs
# PWMS_FLATTENED = os.path.join(PROCESS_DIR, "damo_hg38_PWMS", "flattened")
# # TFBS_ORGANIZED = os.path.join(PROCESS_DIR, "damo_hg38_TFBS")

# # All species IDs for each profile and just those that are human
# PROFILE_MAPPING = os.path.join(PROCESS_DIR, "profile_mapping.txt")
# PROFILE_MAPPING_FILTERED = os.path.join(PROCESS_DIR, "profile_mapping_filtered.txt")

# # Final human specific datasetes
# PWMS_ORGANIZED_FILTER = os.path.join(PROCESS_DIR, "damo_hg38_PWMS_filter")
# TFBS_ORGANIZED_FILTER = os.path.join(PROCESS_DIR, "damo_hg38_TFBS_filter")

# # Biosample maps
# BIOSAMPLE_MAP = os.path.join(PROCESS_DIR, "biosample_map.tsv")

# # Biosamples that failed
# BIOSAMPLE_FAILS = os.path.join(PROCESS_DIR, "biosample_fails.txt")

# # Biosample thresholds
# BIOSAMPLE_THRESHOLDS = os.path.join(PROCESS_DIR, "biosample_thresholds.txt")

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
        EXTRACTED_PWMS,
        EXTRACTED_TFBS,
        # PWMS_ORGANIZED_FILTER,
        # TFBS_ORGANIZED_FILTER,
        # BIOSAMPLE_THRESHOLDS
        #BIOSAMPLE_FAILS,


rule download_unibind:
    message:
        "Downloads Unibind PWMs for TFBs from UniBind database"
    output:
        pwms=PWMS_DOWNLOAD,
        tfbs=TFBS_DOWNLOAD,
    params:
        pwms_url=lambda wc: UNIBIND_URLS["pwms"],
        tfbs_url=lambda wc: UNIBIND_URLS["tfbs"],
    log:
        stdout="workflow/logs/download_unibind.stdout",
        stderr="workflow/logs/download_unibind.stderr",
    conda:
        "../envs/unibind.yaml"
    shell:
        "curl -o {output.pwms} {params.pwms_url} && curl -o {output.tfbs} {params.tfbs_url}"


rule create_mapping:
    message:
        "Creates mapping"
    input:
        rules.download_unibind.output.pwms,
    output:
        BIOSAMPLE_MAP,
    log:
        stdout="workflow/logs/create_mapping.stdout",
        stderr="workflow/logs/create_mapping.stderr",
    conda:
        "../envs/unibind.yaml"
    script:
        "../scripts/mapping.py"


rule flag_filters:
    message:
        "Labels samples to be filtered from analysis"
    input:
        rules.create_mapping.output,
    output:
        mapping=BIOSAMPLE_MAP_FILTER,
        report=BIOSAMPLE_MAP_FILTER + ".report",
    log:
        stdout="workflow/logs/flag_filters.stdout",
        stderr="workflow/logs/flag_filters.stdout",
    conda:
        "../envs/unibind.yaml"
    script:
        "../scripts/filter.py"


rule extract_pwms:
    message:
        "Extracts PWMs from tarball"
    input:
        data=rules.download_unibind.output.pwms,
        mapping=rules.flag_filters.output.mapping,
    output:
        directory(EXTRACTED_PWMS),
    log:
        stdout="workflow/logs/extract_pwms.stdout",
        stderr="workflow/logs/extract_pwms.stdout",
    conda:
        "../envs/unibind.yaml"
    script:
        "../scripts/extract.py"


rule extract_tfbs:
    message:
        "Extracts TFBS from tarball"
    input:
        data=rules.download_unibind.output.tfbs,
        mapping=rules.flag_filters.output.mapping,
    output:
        directory(EXTRACTED_TFBS),
    log:
        stdout="workflow/logs/extract_tfbs.stdout",
        stderr="workflow/logs/extract_tfbs.stdout",
    conda:
        "../envs/unibind.yaml"
    script:
        "../scripts/extract.py"


# rule unpack_unibind_pwms:
#     message:
#         "Unpacks UniBind PWM download"
#     input:
#         rules.download_unibind_pwms.output,
#     output:
#         temp(directory(PWMS_UNPACKED)),
#     params:
#         outdir=lambda w, input: os.path.dirname(input[0]),
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="unpack_unibind_pwms"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="unpack_unibind_pwms"),
#     conda:
#         "../envs/unibind.yaml"
#     shell:
#         "tar -xzf {input} -C {params.outdir}"

# rule flatten_pwms:
#     message:
#         "Flattens PWM dir"
#     input:
#         rules.unpack_unibind_pwms.output,
#     output:
#         directory(PWMS_FLATTENED),
#     params:
#         extension=EXTENSIONS["pwms"],
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="flatten_pwms"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="flatten_pwms"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "../scripts/flatten.py"


# rule organize_pwms:
#     message:
#         "Organizes UniBind damo PWMS by TF and profile"
#     input:
#         rules.unpack_unibind_pwms.output,
#     output:
#         temp(directory(PWMS_ORGANIZED)),
#     params:
#         extension=EXTENSIONS["pwms"],
#         script="../scripts/organize.py",
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="organize_pwms"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="organize_pwms"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "{params.script}"


# rule download_unibind_tfbs:
#     message:
#         "Downloads all TFBSs from UniBind database"
#     output:
#         TFBS_DOWNLOAD,
#     params:
#         url=TFBS_URL,
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="download_unibind_tfbs"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="download_unibind_tfbs"),
#     conda:
#         "../envs/unibind.yaml"
#     shell:
#         "curl -o {output} {params.url}"


# rule unpack_unibind_tfbs:
#     message:
#         "Unpacks UniBind TFBS download"
#     input:
#         rules.download_unibind_tfbs.output,
#     output:
#         temp(directory(TFBS_UNPACKED)),
#     params:
#         outdir=lambda w, input: os.path.dirname(input[0]),
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="unpack_unibind_tfbs"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="unpack_unibind_tfbs"),
#     conda:
#         "../envs/unibind.yaml"
#     shell:
#         "tar -xzf {input} -C {params.outdir}"


# rule organize_tfbs:
#     message:
#         "Organizes UniBind damos by TF and profile"
#     input:
#         rules.unpack_unibind_tfbs.output,
#     output:
#         temp(directory(TFBS_ORGANIZED)),
#     params:
#         extension=EXTENSIONS["tfbs"],
#         script="../scripts/organize.py",
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="organize_tfbs"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="organize_tfbs"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "{params.script}"
# rule profile_mapping:
#     message:
#         "Confirms each profile is homo sapiens."
#     input:
#         rules.organize_pwms.output,
#     output:
#         PROFILE_MAPPING,
#     params:
#         script="../scripts/mapping.py",
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="profile_mapping"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="profile_mapping"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "{params.script}"
# rule filter_tf_targets:
#     message:
#         "Reduces list of profiles to just those that are human and non-redundant."
#     input:
#         rules.profile_mapping.output,
#     output:
#         PROFILE_MAPPING_FILTERED,
#     params:
#         script="../scripts/profiles.py",
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="filter_tf_targets"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="filter_tf_targets"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "{params.script}"
# rule reduce_data_to_filter:
#     message:
#         "Reduces downloaded data to those that are human"
#     input:
#         pwms=rules.organize_pwms.output,
#         tfbs=rules.organize_tfbs.output,
#         profiles=rules.filter_tf_targets.output,
#     output:
#         pwms=directory(PWMS_ORGANIZED_FILTER),
#         tfbs=directory(TFBS_ORGANIZED_FILTER),
#     params:
#         script="../scripts/reduce.py",
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="reduce_data_to_filter"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="reduce_data_to_filter"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "{params.script}"
# rule fetch_thresholds:
#     message:
#         "Creates map of biosamples to recommended score thresholds for scanning."
#     input:
#         rules.reduce_data_to_filter.output.pwms,
#     output:
#         BIOSAMPLE_THRESHOLDS,
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="filter_tf_targets"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="filter_tf_targets"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "../scripts/thresholds.py"
# rule reduct_data_biosamples:
#     message:
#         "Reduces downloaded data to those that have resonable PWMs"
#     input:
#         pwms=rules.organize_pwms.output,
#         tfbs=rules.organize_tfbs.output,
#     output:
#         BIOSAMPLE_FAILS,
#     params:
#         script="../scripts/biosamples.py",
#     log:
#         stdout=expand("workflow/logs/{rule}.stdout", rule="reduct_data_biosamples"),
#         stderr=expand("workflow/logs/{rule}.stderr", rule="reduct_data_biosamples"),
#     conda:
#         "../envs/unibind.yaml"
#     script:
#         "{params.script}"
