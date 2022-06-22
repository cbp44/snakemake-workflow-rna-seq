import pandas as pd

counts = pd.read_table(snakemake.input.tsv, index_col=0)
gene_id_map = pd.read_table(snakemake.input.gene_id_map, index_col=0)

result = pd.concat([id_map,counts],axis=1)

result.to_csv(snakemake.output.tsv, sep="\t")