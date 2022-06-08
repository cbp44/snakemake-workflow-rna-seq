rule RSeQC_GTF2Bed:
    input:
        db=get_ensembl_path("rseqc_annotation.db"),
    output:
        bed="results/qc/rseqc/annotation.bed",
    log:
        "results/logs/rseqc_gtf2bed.log"
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule RSeQC_Junction_Annotation:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.junctionanno.junction.bed",
        "results/qc/rseqc/{sample}-{unit}.junctionanno.junction.Interact.bed"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.log"
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix="results/qc/rseqc/{sample}-{unit}.junctionanno"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule RSeQC_Junction_Saturation:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.log"
    params:
        extra=r"-q 255",
        prefix="results/qc/rseqc/{sample}-{unit}.junctionsat"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule RSeQC_Stat:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}-{unit}.stats.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_stat/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"


rule RSeQC_Infer:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.infer_experiment.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_infer/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule RSeQC_Innerdis:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_innerdis/{sample}-{unit}.log"
    params:
        prefix="results/qc/rseqc/{sample}-{unit}.inner_distance_freq"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule RSeQC_Readdis:
    input:
        bam="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed"
    output:
        "results/qc/rseqc/{sample}-{unit}.readdistribution.txt"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_readdis/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule RSeQC_Readdup:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/qc/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_readdup/{sample}-{unit}.log"
    params:
        prefix="results/qc/rseqc/{sample}-{unit}.readdup"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule RSeQC_ReadGC:
    input:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam"
    output:
        "results/qc/rseqc/{sample}-{unit}.readgc.GC_plot.pdf"
    priority: 1
    log:
        "results/logs/rseqc/rseqc_readgc/{sample}-{unit}.log"
    params:
        prefix="results/qc/rseqc/{sample}-{unit}.readgc"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule MultiQC:
    """
    Create a MultiQC report file using the various rule outputs.
    """
    input:
        expand("results/star/{unit.sample}-{unit.unit}/Aligned.sortedByCoord.out.bam", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.stats.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("results/qc/rseqc/{unit.sample}-{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        expand("results/logs/rseqc/rseqc_junction_annotation/{unit.sample}-{unit.unit}.log", unit=units.itertuples())
    output:
        report("results/qc/multiqc_report.html", "../report/multiqc.rst", category="QC", subcategory="MultiQC")
    log:
        "results/logs/multiqc.log"
    params:
        "" # Optional: extra parameters for multiqc.
    wrapper:
        "0.65.0/bio/multiqc"
