rule Gene_Name_Map:
    """Create a tsv file with a mapping of Ensembl gene ids to gene names.
    """
    input:
        db=get_ensembl_path("gffutils.db"),
    output:
        tsv="results/gene_id_map.tsv",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gene_id_map.py"
