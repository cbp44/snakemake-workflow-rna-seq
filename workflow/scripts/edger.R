log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(edgeR)
# library(biomaRt)

x <- read.delim(snakemake@input[[1]], row.names="gene", check.names=FALSE, stringsAsFactors=FALSE)

# TODO: Get these factors as strings from units.tsv as Snakemake param. Needs to be ordered by 
group <- factor(c(1,1,2,2))
design <- model.matrix(~group)

y <- DGEList(counts=x,group=group)

keep <- filterByExpr(y, design)

y <- y[keep,,keep.lib.sizes=FALSE]

y <- calcNormFactors(y)



y <- estimateDisp(y,design)


et <- exactTest(y)

topTags(et)
