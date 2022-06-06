BEGIN {
    IFS=" "; OFS="\t";
}
{
    gene=$1
    mean_normalized_counts=$2
    log2fc=$3
    stderr=$4
    wald_stat=$5
    wald_test_pvalue=$6
    bh_adjusted_pvalue=$7
    if (bh_adjusted_pvalue<=0.1 && log2fc<0) {
        print gene,mean_normalized_counts,log2fc,bh_adjusted_pvalue
    };
}
