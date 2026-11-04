#!/usr/bin/env python3

import pandas as pd

bedFile=str(snakemake.input.bed_strgm_bdtm)
gtfFile=str(snakemake.input.gtf_strgm_bdtm)
Output=str(snakemake.output.gtf1_out)

Colnames=["Chr","Star","End","gene_id","Length","Strand"]
bedData=pd.read_csv(bedFile, sep="\t", names=Colnames)
bedData["gene_id"]=bedData.gene_id.str.split(";", expand = True)[0]
listId=list(bedData.gene_id)

gtfData=pd.read_csv(gtfFile, sep="\t", header=None)
gtfData["Gene_id"]=gtfData[8].str.split("\"", expand=True)[1]
df=gtfData.query('Gene_id == @listId').drop(["Gene_id"], axis=1)
#print(df)
df.to_csv(Output, sep="\t", header=False, index=False)

