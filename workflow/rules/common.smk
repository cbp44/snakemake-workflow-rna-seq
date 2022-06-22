import os
from snakemake.utils import format as smk_format

def is_human_genome():
    """Returns true if this is configured as a human genome
    """
    return config["ref"]["species"] == "homo_sapiens"

def in_docker_container():
    """Checks if we are running inside a Docker conatiner"""
    with open("/proc/self/cgroup", "r") as procfile:
        for line in procfile:
            fields = line.strip().split("/")
            if "docker" in fields:
                return True
    return False

def get_fastq(wildcards):
    return expand(f"resources/reads/{wildcards.sample}-{wildcards.unit}.{{read_num}}.fastq.gz", read_num=[1,2])
    # return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_ensembl_path(extra_path=""):
    if extra_path != "":
        return "resources/{0}/ensembl/{extra_path}".format(config['ref']['species'], extra_path=extra_path)
    else:
        return "resources/{0}/ensembl".format(config["ref"]["species"])

def get_clinvar_path(extra_path=""):
    if extra_path != "":
        return "resources/{0}/clinvar/{extra_path}".format(config['ref']['species'], extra_path=extra_path)
    else:
        return "resources/{0}/clinvar".format(config["ref"]["species"])

def get_star_index_path():
    ensembl_path = get_ensembl_path()
    return "{ensembl_path}/star_genome".format(ensembl_path=ensembl_path)

def get_star_mane_index_path():
    ensembl_path = get_ensembl_path()
    return "{ensembl_path}/star_genome_mane".format(ensembl_path=ensembl_path)

def get_gtf_annotation_file():
    return get_ensembl_path("genome.gtf.gz")
    # return smk_format("resources/{config['ref']['organism']}/ensembl/genome.gtf.gz")

def is_single_end(sample, unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]


def get_fq(wildcards, read_num=1, host_path=True):
    if host_path:
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1_full" if read_num == 1 else "fq2_full"]].dropna()
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return f"resources/reads/{wildcards.sample}-{wildcards.unit}.{read_num}.fastq.gz"
    else:
        # yes trimming, use trimmed data
        return expand("results/trimmed/{sample}-{unit}.{read_num}.fastq.gz",
                      read_num=read_num, **wildcards)


def get_fq1(wildcards):
    return get_fq(wildcards, read_num=1, host_path=False)


def get_fq2(wildcards):
    return get_fq(wildcards, read_num=2, host_path=False)

def get_fq1_host(wildcards, host_path=True):
    return get_fq(wildcards, read_num=1, host_path=True)

def get_fq2_host(wildcards, host_path=True):
    return get_fq(wildcards, read_num=2, host_path=True)

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6
