localrules: download_files_local

def find_unit_by_fastq_file(filename, read_col="fq1"):
    if filename.endswith("1"):
        fq_col = "fq1"
        fq_full_col = "fq1_full"
    elif filename.endswith("2"):
        fq_col = "fq2"
        fq_full_col = "fq2_full"

    units_reindexed = units.set_index(fq_col)

    found_files = units_reindexed.loc[f"resources/reads/{filename}.fq.gz", [fq_full_col]]

    return found_files

rule download_files_local:
    input:
        lambda wc: find_unit_by_fastq_file(wc.filename)
    output:
        "resources/reads/{filename}.fq.gz"
    shell:
        "cp -av --dereference {input} {output}"