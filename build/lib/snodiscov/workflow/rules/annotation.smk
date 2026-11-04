# =============================================================================
#  Build the PIC (potential interesting candidates) annotation and the two
#  extended annotations used for quantification:
#    annotFile1 = reference GTF + all PICs           -> CoCo on discovery data
#    annotFile2 = reference GTF + expressed PICs     -> STAR/CoCo on validation data
# =============================================================================
PA = "work/pics_annotation"
ANN = "work/annotation"
GOI = "results/gene_of_interest"   # main output 1: expressed candidate transcripts


localrules:
    concat_strgm_bdtm,
    annotFile1,
    annotFile2,
    concat_pics_bdg_parsing,


rule concat_strgm_bdtm:
    input:
        strgm=rules.annotFile1_strgM.output[0],
        bdtm=rules.annotFile1_bdtM.output[0],
    output:
        PA + "/pics_annotation_file1_project1_strgm_bdtm.gtf",
    conda:
        ENVS + "/coreutils.yaml"
    shell:
        "cat {input.strgm} {input.bdtm} > {output}"


rule merge_strgm_bdtm:
    input:
        rules.concat_strgm_bdtm.output[0],
    output:
        bed=PA + "/pics_annotation_file1_project1_strgm_bdtm.bed",
    conda:
        ENVS + "/coco.yaml"
    log:
        "logs/discovery/merge_strgm_bdtm.log",
    shell:
        """
        awk '$3=="gene" {{OFS="\\t"; print ($1,$4,$5,$10,$5-$4,$7)}}' {input} \
            | bedtools sort -i stdin \
            | bedtools merge -s -i stdin -c 4 -o collapse -delim " " > {output.bed} 2> {log}
        """


rule picsAnnotFile1:
    input:
        # kept as a one-element list: same object type as in the original rule
        bed_strgm_bdtm=[rules.merge_strgm_bdtm.output.bed],
        gtf_strgm_bdtm=rules.concat_strgm_bdtm.output[0],
    output:
        gtf1_out=PA + "/pics_annotation_file1_project1_tmp.gtf",
    log:
        "logs/discovery/picsAnnotFile1.log",
    conda:
        ENVS + "/python.yaml"
    script:
        "../scripts/getGtf_strgm_bdtm.py"


rule parsPicsAnnotFile1:
    """Remove double quotes in the GTF frame field."""
    input:
        picgtf=rules.picsAnnotFile1.output.gtf1_out,
    output:
        parsGtfPics=PA + "/pics_annotation_file1_project1.gtf",
    params:
        script=SCRIPTS + "/parsAnnotFile1.sh",
    log:
        "logs/discovery/parsPicsAnnotFile1.log",
    conda:
        ENVS + "/coreutils.yaml"
    shell:
        "bash {params.script} {input.picgtf} {output.parsGtfPics} &> {log}"


rule annotFile1:
    input:
        gtfRef=REF_GTF,
        gtfPics=rules.parsPicsAnnotFile1.output.parsGtfPics,
    output:
        ANN + "/hg38_Ensembl_V109_Scottlab_2024_snoRNAs_tRNA_pics.gtf",
    conda:
        ENVS + "/coreutils.yaml"
    shell:
        "cat {input.gtfRef} {input.gtfPics} > {output}"


rule coco_ca:
    """CoCo correct_annotation (writes <name>.correct_annotation.gtf next to input)."""
    input:
        gtf=ANN + "/{name}.gtf",
        coco=COCO_DONE,
    output:
        ANN + "/{name}.correct_annotation.gtf",
    wildcard_constraints:
        name=r"[^/]+",
    params:
        coco_bin=COCO_BIN,
    resources:
        mem_mb=30000,
        runtime=180,
    log:
        "logs/annotation/{name}_coco_ca.log",
    conda:
        ENVS + "/coco.yaml"
    shell:
        "export PATH={params.coco_bin}:$PATH; coco ca {input.gtf} &> {log}"


rule coco_parsing:
    """Keep PICs detected by CoCo in the discovery samples (GTF + BED)."""
    input:
        tab=expand("work/{d}/coco_cc/{s}.tsv", d=D, s=DISC_SAMPLES),
        gtf=rules.parsPicsAnnotFile1.output.parsGtfPics,
    output:
        gtf_pics=GOI + "/gene_of_interest.gtf",
        bed_pics=GOI + "/gene_of_interest.bed",
    params:
        folder="work/" + D + "/coco_cc/",
        script=SCRIPTS + "/cocoParsingCount.py",
    conda:
        ENVS + "/python.yaml"
    log:
        "logs/discovery/coco_parsing.log",
    shell:
        """
        python {params.script} --gtffile {input.gtf} --folder {params.folder} \
            --run_mode gtffile > {output.gtf_pics} 2> {log}
        python {params.script} --gtffile {input.gtf} --folder {params.folder} \
            --run_mode bedfile > {output.bed_pics} 2>> {log}
        """


rule annotFile2:
    input:
        gtf_hg38=REF_GTF,
        gtf_pics=rules.coco_parsing.output.gtf_pics,
    output:
        ANN + "/hg38_Ensembl_V109_Scottlab_2024_snoRNAs_pics_dep.gtf",
    conda:
        ENVS + "/coreutils.yaml"
    shell:
        "cat {input.gtf_hg38} {input.gtf_pics} > {output}"


# ------------------------------------------------- coverage around PICs ----
rule pics_bdg_intersect:
    input:
        bdg=lambda wc: "work/" + D + "/coco_cb/" + str(wc.sample) + "_sorted.bedgraph",
        bed=rules.coco_parsing.output.bed_pics,
    output:
        "work/discovery/bedgraph_filter/intersect_out/{sample}_sorted_bedgraph.bed",
    params:
        script=SCRIPTS + "/bedgraph_filtering.sh",
    conda:
        ENVS + "/coco.yaml"
    log:
        "logs/discovery/bdg_intersect/{sample}.log",
    shell:
        "bash {params.script} {input.bdg} {input.bed} {output} &> {log}"


rule pics_bdg_parsing:
    input:
        rules.pics_bdg_intersect.output[0],
    output:
        "work/discovery/bedgraph_filter/parsing_out/{sample}_sorted_bedgraph_parsing.tsv",
    params:
        script=SCRIPTS + "/bedgraph_filtering.py",
    threads: 8
    resources:
        mem_mb=80000,
        runtime=480,
    conda:
        ENVS + "/python.yaml"
    log:
        "logs/discovery/bdg_parsing/{sample}.log",
    shell:
        "python {params.script} --bdg_data {input} --output_gene {output} &> {log}"


rule concat_pics_bdg_parsing:
    input:
        expand(
            "work/discovery/bedgraph_filter/parsing_out/{s}_sorted_bedgraph_parsing.tsv",
            s=sorted(DISC_SAMPLES),
        ),
    output:
        "work/discovery/bedgraph_filter/pics_bedgraph_filter_project1.tsv",
    conda:
        ENVS + "/coreutils.yaml"
    shell:
        "cat {input} > {output}"
