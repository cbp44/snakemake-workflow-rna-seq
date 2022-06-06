import os
from snakemake.utils import format as smk_format

def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    
def get_ensembl_path(extra_path=""):
    if extra_path != "":
        return "resources/{0}/ensembl/{extra_path}".format(config['ref']['organism'], extra_path=extra_path)
    else:
        return "resources/{0}/ensembl".format(config["ref"]["organism"])

def get_star_index_path():
    ensembl_path = get_ensembl_path()
    return "{ensembl_path}/star_index".format(ensembl_path=ensembl_path)

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


def get_fq(wildcards, read_num=1):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1" if read_num == 1 else "fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        return expand("results/trimmed/{sample}-{unit}.{read_num}.fastq.gz",
                      read_num=read_num, **wildcards)


def get_fq1(wildcards):
    return get_fq(wildcards, read_num=1)


def get_fq2(wildcards):
    return get_fq(wildcards, read_num=2)



def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6
