from collections import defaultdict
import csv 

def get_genes_from_exon_bed(bed_filename, fieldnames=None):
    # Each gene has a list containing the exon rows from the BED file
    genes = defaultdict(lambda : [])
    with open(bed_filename, 'r') as bed_file:
        if fieldnames:
            reader = csv.DictReader(bed_file, delimiter='\t', fieldnames=fieldnames)
        else:
            reader = csv.DictReader(bed_file, delimiter='\t')
        
        for row in reader:
            gene_hgnc = row.get('gene_hgnc', None)
            if gene_hgnc:
                # genes[gene_hgnc].append(row)
                genes[row['ensembl']].append(row)

    return genes


def get_gene_list(filename):
    with open(filename, 'r') as gene_list_file:
        reader = csv.reader(gene_list_file)
        gene_list = [row[0] for row in reader if len(row) > 0]

    return gene_list


fieldnames=["chr","start","end","rank","ensembl","refseq","gene_hgnc"]


def write_gene_exons(bed_filename, gene_exons):
    with open(bed_filename, 'w') as bed_file:
        writer = csv.DictWriter(bed_file, delimiter='\t', fieldnames=fieldnames)
       
        for exon in gene_exons:
            writer.writerow(exon)


# Get the exons of every gene as a list in a dict where the key is the refseq accession number
all_genes = get_genes_from_exon_bed(snakemake.input.cds_regions)

# Pass in the field names here because the input file was created by bedtools which removes the header
target_genes = get_gene_list(snakemake.input.diffexp_genes)

gene_exons_to_keep = []
# gene_exons_to_discard = []

# print(all_genes)
# print(target_genes)

for gene_name in target_genes:
    # print(gene_name)
    if gene_name in all_genes.keys():
        gene_exons_to_keep.extend(all_genes[gene_name])

write_gene_exons(snakemake.output.bed, gene_exons_to_keep)
