rule Count_Matrix:
    input:
        expand("results/star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "results/counts/all.tsv"
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units)
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

use rule Count_Matrix as Count_Matrix_MANE with:
    input:
        expand("results/star_mane/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "results/counts_mane/all.tsv"


rule Map_Gene_IDs:
    input:
        tsv="results/{path}/{filename}.tsv",
        gene_id_map="results/gene_id_map.tsv",
    output:
        tsv="results/{path}/{filename}.with_gene_name.tsv",
    conda:
        "../envs/pandas.yaml"
    wildcard_constraints:
        counts_path="|".join(["counts","counts_mane","diffexp","diffexp_mane","deseq2","deseq2_mane"])
    script:
        "../scripts/map_gene_ids.py"

rule DESeq2_Init:
    input:
        counts="results/counts/all.tsv"
    output:
        "results/deseq2/all.rds",
        "results/deseq2/normcounts.tsv",
    params:
        samples=config["samples"],
        design=config["diffexp"]["model"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "results/logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

use rule DESeq2_Init as DESeq2_Init_MANE with:
    input:
        counts="results/counts_mane/all.tsv"
    output:
        "results/deseq2_mane/all.rds",
        "results/deseq2_mane/normcounts.tsv",
    log:
        "results/logs/deseq2_mane/init.log"

rule PCA:
    input:
        "results/deseq2/all.rds"
    output:
        report("results/pca.svg", "../report/pca.rst", category="QC", subcategory="PCA")
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "results/logs/pca.log"
    script:
        "../scripts/plot-pca.R"


rule DESeq2:
    input:
        "results/deseq2/all.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", 
            "../report/diffexp.rst", category="Differential Expression", subcategory="{contrast}"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst", category="Differential Expression", subcategory="{contrast}"),
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "results/logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"


use rule DESeq2 as DESeq2_MANE with:
    input:
        "results/deseq2_mane/all.rds"
    output:
        table=report("results/diffexp_mane/{contrast}.diffexp.tsv", 
            "../report/diffexp.rst", category="Differential Expression (MANE)", subcategory="{contrast}"),
        ma_plot=report("results/diffexp_mane/{contrast}.ma-plot.svg", "../report/ma.rst", category="Differential Expression (MANE)", subcategory="{contrast}"),
    log:
        "results/logs/deseq2_mane/{contrast}.diffexp.log"


wildcard_constraints:
    diffexp="|".join(["diffexp","diffexp_mane"])

rule Top_Upregulated_Genes:
    input:
        "results/{diffexp}/{contrast}.diffexp.tsv",
    output:
        "results/{diffexp}/{contrast}.diffexp.top_upregulated.tsv"
    conda: "../envs/awk.yml"
    shell:
        """
        (head -n1 {input}; \
        tail -n+2 {input} | awk 'BEGIN {{OFS="\\t"}} {{if ($7<=0.1 && $3>0) {{ \
                print $1,$2,$3,$4,$5,$6,$7}} \
            }}' | sort -k5,5nr -k3,3n) \
        > {output}
        """

rule Top_Downregulated_Genes:
    input:
        "results/{diffexp}/{contrast}.diffexp.tsv",
    output:
        "results/{diffexp}/{contrast}.diffexp.top_downregulated.tsv"
    conda: "../envs/awk.yml"
    shell:
        """
        (head -n1 {input}; \
        tail -n+2 {input} | awk 'BEGIN {{OFS="\\t"}} {{if ($7<=0.1 && $3<0) {{ \
                print $1,$2,$3,$4,$5,$6,$7}} \
            }}' | sort -k5,5nr -k3,3n) \
        > {output}
        """


# use rule Top_Upregulated_Genes as Top_Upregulated_MANE_Genes with:
#     input:
#         "results/diffexp_mane/{contrast}.diffexp.tsv",
#     output:
#         "results/diffexp_mane/{contrast}.diffexp.top_upregulated.tsv"


# use rule Top_Downregulated_Genes as Top_Downregulated_MANE_Genes with:
#     input:
#         "results/diffexp_mane/{contrast}.diffexp.tsv",
#     output:
#         "results/diffexp_mane/{contrast}.diffexp.top_downregulated.tsv"
