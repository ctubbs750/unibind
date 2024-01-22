from os import listdir, path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

PWMS_URL = config["urls"]["pwms"]
TFBS_URL = config["urls"]["tfbs"]
OUTP_DIR = config["install_dir"]

# ------------- #
# I/O           #
# ------------- #

# Raw PWM and TFBS download
PWMS_DOWNLOAD = os.path.join(OUTP_DIR, "damo_hg38_PWMs.tar.gz")
TFBS_DOWNLOAD = os.path.join(OUTP_DIR, "damo_hg38_TFBS_per_TF.tar.gz")

# Unpacked PWMs
PWMS_UNPACKED = os.path.join(OUTP_DIR, "damo_hg38_PWMs")
TFBS_UNPACKED = os.path.join(OUTP_DIR, "damo_hg38_TFBS_per_TF")

# Orgnaized PWMs
PWMS_ORGANIZED = os.path.join(OUTP_DIR, "damo_hg38_PWMS")
TFBS_ORGANIZED = os.path.join(OUTP_DIR, "damo_hg38_TFBS")


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
    default_target: True


rule download_unibind_pwms:
    message:
        """
        Downloads all PWMS from UniBind database
        """
    output:
        temp(PWMS_DOWNLOAD),
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
        directory(PWMS_ORGANIZED),
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
        temp(TFBS_DOWNLOAD),
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
        outdir=OUTP_DIR,
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
        directory(TFBS_ORGANIZED),
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
