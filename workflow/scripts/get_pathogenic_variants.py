## Author: Christopher B. Preusch
## Date: Nov 8, 2021

from collections import defaultdict
import csv 
import vcf

# Initialize a defaultdict for genes to hold MC with genes
genes = defaultdict(lambda : defaultdict(lambda : 0))

# This will hold all possible molecular consequence values
# to be used for printing the final table
possible_mcs = []

def process_vcf(vcf_filename):
    with open(vcf_filename, 'r') as vcf_file:
        vcf_reader = vcf.Reader(vcf_file)
        for record in vcf_reader:
            clnsig = record.INFO.get("CLNSIG", [])

            # Only process record if it is pathogenic
            if "Pathogenic" in clnsig or "Pathogenic/Likely_pathogenic" in clnsig or "Likely_pathogenic" in clnsig:
        
                molecular_consequences = record.INFO.get("MC", [])
                molecular_consequences = [x.split("|")[1] for x in molecular_consequences]
                possible_mcs.extend(molecular_consequences)
            
                # Only process records with a gene attached
                gene_info = record.INFO.get("GENEINFO", None)
                if gene_info:
                    gene = gene_info.split(":")[0]
                    for mc in molecular_consequences:
                        genes[gene][mc] += 1
                    
            
def write_tsv(mcs, genes_dict, tsv_filename):
    # The header of the tsv file
    header = ["gene_name"]
    header.extend(mcs)

    with open(tsv_filename, 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(header)
        for gene, mc_dict in genes_dict.items():
            gene_line = [gene]
            for mc in mcs:
                mc_count = mc_dict.get(mc, 0)
                gene_line.append(mc_count)
                    
            writer.writerow(gene_line)



process_vcf(snakemake.input[0])

# Create a unique set of possible molecular consequences

possible_mcs = list(set(possible_mcs))
possible_mcs.sort()

write_tsv(possible_mcs, genes, snakemake.output[0])
