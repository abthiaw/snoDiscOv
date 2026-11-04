#!/usr/bin/env python3

import picsParsing

## 4. Prsing NOP58 and DKC1 depletion results from DESeq2

picsParsing.deseq2Parsing(
    deseq2tab=snakemake.input.deseq2file,
    snoGtf=snakemake.input.snogtf,
    tab_out=snakemake.output.DepOut,
    testType=snakemake.params.depType
)