log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tidyverse")
library("DESeq2")
#library(EnhancedVolcano)
#library(sjPlot)
#library(readxl)
#library(reshape2)


dataCount=read_tsv(snakemake@input[["Count"]])%>%
  column_to_rownames(var="gene_id")
dataCount=dataCount[,order(names(dataCount))]
i=c(1:ncol(dataCount))
dataCount[ , i] <- apply(dataCount[ , i], 2, function(x) as.integer(as.character(x)))
metadata=read_tsv(snakemake@input[["condition"]])%>%
    arrange(samples)
# Align count columns on the design table by name (the original relied on
# order() and arrange() sorting identically, which depends on the locale)
stopifnot(setequal(colnames(dataCount), metadata$samples))
dataCount=dataCount[, metadata$samples]
metadata$features=as.factor(metadata$features)
metadata$samples=as.factor(metadata$samples)

#Run DESeq2
rbp_dds=DESeqDataSetFromMatrix(countData=dataCount, colData=metadata, design=~features)
rbp_dds$features = relevel(rbp_dds$features, "control")
dds=DESeq(rbp_dds)
res=results(dds)
dfres=as.data.frame(res)%>%rownames_to_column() #Saving  res as a dataframe
dfres=dfres[complete.cases(dfres), ] # Remove rows with NAN values

## Define 2 genes group list (annotated (an_gene) or non annotated (nan_gene)) 
nan_gene_data=dfres%>%filter(grepl("STRG", rowname)) #Non annotated gene
colnames(nan_gene_data)=c('Gene_Id','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
write.table(nan_gene_data, file=snakemake@output[["deseq2_out"]], sep='\t', row.names =FALSE, col.names=TRUE)
