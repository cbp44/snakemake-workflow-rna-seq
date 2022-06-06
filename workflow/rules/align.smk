rule align:
    input:
        fq1 = get_fq1,
    	fq2 = get_fq2
    output:
        # see STAR manual for additional output files
        "results/star/{sample}-{unit}/Aligned.out.bam",
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/{sample}-{unit}/Log.final.out",
        "results/star/{sample}-{unit}/Log.out",
        "results/star/{sample}-{unit}/Log.progress.out",
        "results/star/{sample}-{unit}/ReadsPerGene.out.tab",
        "results/star/{sample}-{unit}/SJ.out.tab"
    shadow: "shallow"
    log:
        "results/logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index=get_star_index_path(),
        # optional parameters
        extra="--outSAMtype BAM Unsorted SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {0} {1}".format(get_gtf_annotation_file(), config['params']['star']),
    resources:
        mem_mb=128000
    threads: 24
    wrapper:
        "0.65.0/bio/star/align"

rule samtools_index:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam.bai"
    params:
        "" # optional params string
    resources:
        mem_mb=6000
    wrapper:
        "0.65.0/bio/samtools/index"

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
