wildcard_constraints:
    regulation_direction="|".join(["upregulated","downregulated"])


rule Get_DE_Gene_List:
    input:
        "results/diffexp/{contrast}.diffexp.top_{regulation_direction}.tsv"
    output:
        temp("results/variants/{contrast}.{regulation_direction}_genes.tsv"),
    conda: "../envs/awk.yml"
    shell:
        "cut -f 1 {input} | sed 's/\"//g' | tail -n+2 > {output}"


rule Get_Gene_CDS_Regions:
    input:
        diffexp_genes="results/variants/{contrast}.{regulation_direction}_genes.tsv",
        all_exons=get_ensembl_path("cds_regions.tsv"),
    output:
        exon_bed="results/variants/{contrast}.{regulation_direction}_gene_exons.bed"
    conda: "../envs/awk.yml"
    script:
        "../scripts/get_diffexp_exons.py"


# rule Fix_Chromosome_Prefixes:
#     """
#     Remove the 'chr' chromosome prefixes in the min max coords BED file to
#     match the ones used in the ClinVar VCF file.
#     """
#     input:
#         "results/variants/{contrast}.{regulation_direction}_gene_exons.bed"
#     output:
#         temp("results/variants/{contrast}.{regulation_direction}_gene_exons.no_prefixes.bed")
#     version: "0.1"
#     conda: "../envs/awk.yml"
#     shell:
#         "sed s/^chr//g {input[0]} > {output[0]}"

rule Filter_Bad_Variants:
    input:
        get_ensembl_path("clinically_associated_variants.vcf.gz")
    output:
        get_ensembl_path("clinically_associated_variants.filtered.vcf.gz")
    conda: "../envs/bcftools.yaml"
    shell:
        "bcftools view -V other {input} | gzip - > {output}"

rule Variants_in_CDS_Regions:
    """
    Get all variants in the input.clinvar_vcf file that are within input.target_region_bed.
    Output is VCF subset of the original VCF.
    """
    input:
        target_region_bed="results/variants/{contrast}.{regulation_direction}_gene_exons.bed",
        vcf=get_ensembl_path("clinically_associated_variants.vcf.gz"),
        # clinvar_vcf=get_ensembl_path("clinvar_variants.vcf.gz"),
    output:
        "results/variants/{contrast}.clinically_assoc_variants_in_{regulation_direction}_gene_exons.vcf"
    version: "0.1"
    conda: "../envs/bedtools.yml"
    shell:
        "zcat {input.vcf} | bedtools intersect -a stdin -b {input.target_region_bed} -header -wa -u > {output}"


rule Extract_Pathogenic_Variants:
    """
    Gets the variants in a VCF file classified with keys "Pathogenic", 
    "Likely_pathogenic" or "Pathogenic/Likely_pathogenic"
    """
    input:
        "results/variants/{contrast}.clinically_assoc_variants_in_{regulation_direction}_gene_exons.vcf"
    output:
        "results/variants/{contrast}.{regulation_direction}_gene_pathogenic_variant_summary.tsv"
    version: "0.1"
    conda: "../envs/pyvcf.yml"
    script: "../scripts/get_pathogenic_variants.py"
