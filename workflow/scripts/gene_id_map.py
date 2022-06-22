import gffutils

db = gffutils.FeatureDB(snakemake.input.db)

id_name_map = dict()

with open(snakemake.output.tsv, 'w') as outfileobj:
    outfileobj.write('{}\n'.format('\t'.join(['gene_id','gene_name'])))
    for gene in db.features_of_type('gene'):
        if "gene_id" in gene.attributes: 
            if "gene_name" in gene.attributes:
                id_name_map[gene.attributes["gene_id"][0]] = gene.attributes["gene_name"][0]
            else:
                id_name_map[gene.attributes["gene_id"][0]] = gene.attributes["gene_id"][0] 
    for gene_id,gene_name in id_name_map.items():
        outfileobj.write('{}\n'.format('\t'.join([gene_id,gene_name])))