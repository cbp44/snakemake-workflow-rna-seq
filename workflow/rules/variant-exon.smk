wildcard_constraints:
    regulation_direction="|".join(["upregulated","downregulated"])


rule Get_DE_Gene_List:
    input:
        "results/diffexp_mane/{contrast}.diffexp.top_{regulation_direction}.tsv"
    output:
        temp("results/variants/{contrast}.{regulation_direction}_genes.tsv"),
    conda: "../envs/awk.yml"
    shell:
        "cut -f 1 {input} | sed 's/\"//g' | tail -n+2 > {output}"


rule Get_Gene_CDS_Regions:
    input:
        diffexp_genes="results/variants/{contrast}.{regulation_direction}_genes.tsv",
        cds_regions=get_ensembl_path("mane-cds_regions.tsv"),
    output:
        bed="results/variants/{contrast}.{regulation_direction}_gene_cds_regions.bed"
    conda: "../envs/awk.yml"
    script:
        "../scripts/get_diffexp_cds.py"


rule Variants_in_CDS_Regions:
    """
    Get all variants in the input.clinvar_vcf file that are within input.target_region_bed.
    Output is VCF subset of the original VCF.
    """
    input:
        target_region_bed="results/variants/{contrast}.{regulation_direction}_gene_cds_regions.bed",
        # vcf=get_ensembl_path("clinically_associated_variants.filtered.vcf.gz"),
        vcf=get_clinvar_path("clinvar.vcf.gz"),
    output:
        "results/variants/{contrast}.clinically_assoc_variants_in_{regulation_direction}_gene_cds.vcf"
    version: "0.1"
    conda: "../envs/bedtools.yml"
    shell:
        "zcat {input.vcf} | bedtools intersect -a stdin -b {input.target_region_bed} -header -wa -u > {output}"


rule Summarize_Pathogenic_Variant_MCs:
    """
    Gets the variants in a VCF file classified with keys "Pathogenic", 
    "Likely_pathogenic" or "Pathogenic/Likely_pathogenic"
    """
    input:
        "results/variants/{contrast}.clinically_assoc_variants_in_{regulation_direction}_gene_cds.vcf"
    output:
        "results/variants/{contrast}.{regulation_direction}_gene_pathogenic_variant_summary.tsv"
    version: "0.1"
    conda: "../envs/pyvcf.yml"
    script: "../scripts/get_pathogenic_variants.py"
