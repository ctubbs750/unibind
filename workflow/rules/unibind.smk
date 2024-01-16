from snakemake.utils import min_version


# Configuration
configfile: "config/config.yaml"


# Parameters TODO: think about how this working when doing the same in scan...
PWMS_URL = config["urls"]["pwms"]
TFBS_URL = config["urls"]["tfbs"]

# Settings
min_version("7.32.4")


rule all:
    input:
        "resources/data/unibind/damo_hg38_PWMS",
        "resources/data/unibind/damo_hg38_TFBS",
    default_target: True


rule download_unibind_pwms:
    message:
        """
        Downloads all PWMS from UniBind database
        """
    output:
        temp("resources/data/unibind/damo_hg38_PWMs.tar.gz"),
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
        temp(directory("resources/data/unibind/unibind/damo_hg38_PWMs_unpacked")),
    log:
        stdout="workflow/logs/unpack_unibind_pwms.stdout",
        stderr="workflow/logs/unpack_unibind_pwms.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        """
        mkdir -p {output} && tar -xzf {input} -C {output}
        """


rule flatten_unibind_pwms:
    message:
        """
        Flattens unpacked tarbell. Not great - for the moment just check for creation of first in OP.
        Also removes all the empty dirs.
        """
    input:
        rules.unpack_unibind_pwms.output,
    output:
        directory("resources/data/unibind/damo_hg38_PWMS"),
    log:
        stdout="workflow/logs/flatten_unibind_pwms.stdout",
        stderr="workflow/logs/flatten_unibind_pwms.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    shell:
        """
        mkdir {output} &&
        find {input} -mindepth 2 -type f -exec mv -t {output} -i '{{}}' + &&
        find {input} -type d -empty -delete
        """


rule download_unibind_tfbs:
    message:
        """
        Downloads all TFBSs from UniBind database
        """
    output:
        temp("resources/data/unibind/damo_hg38_TFBS_per_TF.tar.gz"),
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
        temp(directory("resources/data/unibind/damo_hg38_TFBS_per_TF")),
    params:
        outdir="resources/data/unibind"
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


rule organize_files:
    message:
        """
        Orgnzies UniBind damos by TF and profile
        """
    input:
        rules.unpack_unibind_tfbs.output,
    output:
        directory("resources/data/unibind/damo_hg38_TFBS"),
    log:
        stdout="workflow/logs/organize_files.stdout",
        stderr="workflow/logs/organize_files.stderr",
    conda:
        "../envs/unibind.yaml"
    threads: 1
    script:
        "../scripts/organize.py"