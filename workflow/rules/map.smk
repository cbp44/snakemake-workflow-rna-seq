rule Copy_Reads:
    input:
        fq1 = get_fq1_host,
        fq2 = get_fq2_host,
    output:
        fq1 = "resources/reads/{sample}-{unit}.1.fastq.gz",
        fq2 = "resources/reads/{sample}-{unit}.2.fastq.gz"
    threads: 6
    shell:
        "cp -v {input.fq1} {output.fq1}; "
        "cp -v {input.fq2} {output.fq2}; "
        "chmod 0644 {output.fq1} {output.fq2}"


rule Map_Reads:
    input:
        fq1 = get_fq1,
        fq2 = get_fq2,
    output:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        log="results/star/{sample}-{unit}/Log.out",
        log_final="results/star/{sample}-{unit}/Log.final.out",
        log_progress="results/star/{sample}-{unit}/Log.progress.out",
        reads_per_gene="results/star/{sample}-{unit}/ReadsPerGene.out.tab",
        sj="results/star/{sample}-{unit}/SJ.out.tab",
    log:
        "results/logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        idx=get_star_index_path(),
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"
    threads: 24
    resources:
        tmpdir=lambda wc: "/bigtmp" if in_docker_container() else "/tmp"
    wrapper:
        "v1.5.0/bio/star/align"

use rule Map_Reads as Map_Reads_MANE with:
    output:
        bam="results/star_mane/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        log="results/star_mane/{sample}-{unit}/Log.out",
        log_final="results/star_mane/{sample}-{unit}/Log.final.out",
        log_progress="results/star_mane/{sample}-{unit}/Log.progress.out",
        reads_per_gene="results/star_mane/{sample}-{unit}/ReadsPerGene.out.tab",
        sj="results/star_mane/{sample}-{unit}/SJ.out.tab",
    log:
        "results/logs/star_mane/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        idx=get_star_mane_index_path(),
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts"


rule Samtools_Index:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai"
    params:
        "" # optional params string
    resources:
        mem_mb=6000
    wrapper:
        "v1.5.0/bio/samtools/index"

rule signal_track:
    """
    Create a normalized signal track for display on the UCSC genome browser.
    """
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai",
    output:
        "results/tracks/signal/{sample}-{unit}.bw"
    log:
        "results/logs/tracks/signal/{sample}-{unit}.log"
    threads: 12
    resources:
        mem_mb=10000
    conda:
        "../envs/deeptools.yaml"
    shell: "bamCoverage --normalizeUsing RPKM --verbose -b {input[0]} -o {output[0]} -p {threads} &> {log}"
