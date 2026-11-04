Project scripts called by the workflow (python/bash are invoked explicitly,
no chmod needed). `picsParsing.py` is a module imported by DepParsing.py
(Snakemake adds this directory to sys.path for `script:` rules).

| script | rule |
|---|---|
| stringtie_geneId.py | stringtie_geneId |
| getAnnotFileProject1.py | target_trans_bdtM / target_trans_strgM |
| buildAnnotFileBdtm.py / buildAnnotFileStrgm.py | annotFile1_bdtM / annotFile1_strgM |
| getGtf_strgm_bdtm.py, parsAnnotFile1.sh | picsAnnotFile1, parsPicsAnnotFile1 |
| cocoParsingCount.py, cocoCount_merged.py | coco_parsing, coco_merge |
| bedgraph_filtering.sh / .py | pics_bdg_intersect / pics_bdg_parsing (coverage_tracks only) |
| deseq2_analysis_smk.R | deseq2 |
| DepParsing.py (+ picsParsing.py), snu13_snocand.py | dep_parsing, snu13_parsing |
