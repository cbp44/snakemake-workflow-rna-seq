BEGIN {
    IFS="\t"; OFS="\t";
}
{
    n=split($10, exon_starts, ",");
    split($11, exon_ends, ",");
    chromosome=$3
    split($2, gene_accession_nr, ".");
    common_gene_name=$13
    strand=$4

    for (i=1; i<n; i++) {
        print chromosome,exon_starts[i],exon_ends[i],"-","-",gene_accession_nr[1],common_gene_name
    };
}
