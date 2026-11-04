# =============================================================================
#  Discovery of candidate transcripts (PICs) from the discovery dataset
#  Two complementary strategies, as in the original pipeline:
#    bdtm  : per-sample StringTie -> gffcompare -> bedtools merge across samples
#    strgm : StringTie --merge across samples -> gffcompare
# =============================================================================
rule stringtie:
    input:
        gtf=REF_GTF,
        bam=lambda wc: bam(D, wc.sample),
    output:
        "work/discovery/stringtie/{sample}_strg.gtf",
    params:
        args=config.get("stringtie", "--fr -m 50"),
    threads: 32
    resources:
        mem_mb=100000,
        runtime=360,
    log:
        "logs/discovery/stringtie/{sample}.log",
    conda:
        ENVS + "/stringtie.yaml"
    shell:
        "stringtie {params.args} -p {threads} -G {input.gtf} -o {output} {input.bam} &> {log}"


# ------------------------------------------------------------------ bdtm ----
rule stringtie_geneId:
    """Make StringTie gene_ids unique across samples. The original used
    $SLURM_JOB_ID as suffix (non-reproducible); now deterministic."""
    input:
        "work/discovery/stringtie/{sample}_strg.gtf",
    output:
        "work/discovery/bdtm/{sample}_stringtie.gtf",
    params:
        suffix=_geneid_suffix,
        script=SCRIPTS + "/stringtie_geneId.py",
    log:
        "logs/discovery/bdtm/{sample}_geneId.log",
    conda:
        ENVS + "/python.yaml"
    shell:
        "python {params.script} --suffixe {params.suffix} --input {input} > {output} 2> {log}"


rule gffcomp_bdtM:
    input:
        gtf_ref=REF_GTF,
        gtf_strg="work/discovery/bdtm/{sample}_stringtie.gtf",
    output:
        multiext("work/discovery/gffcomp_bdtm/{sample}", ".annotated.gtf", ".loci", ".stats", ".tracking"),
        # gffcompare writes .tmap/.refmap next to the *input* file
        tmap="work/discovery/bdtm/{sample}.{sample}_stringtie.gtf.tmap",
    params:
        prefix="work/discovery/gffcomp_bdtm/{sample}",
    log:
        "logs/discovery/gffcomp_bdtm/{sample}.log",
    conda:
        ENVS + "/gffcompare.yaml"
    shell:
        "gffcompare -R -r {input.gtf_ref} -o {params.prefix} {input.gtf_strg} &> {log}"


rule target_trans_bdtM:
    input:
        tmap="work/discovery/bdtm/{sample}.{sample}_stringtie.gtf.tmap",
        gtf_strg="work/discovery/bdtm/{sample}_stringtie.gtf",
    output:
        gtf="work/discovery/pics/{sample}_bdtm.gtf",
        tab="work/discovery/pics/{sample}_bdtm.tsv",
    params:
        script=SCRIPTS + "/getAnnotFileProject1.py",
    log:
        "logs/discovery/pics_bdtm/{sample}.log",
    conda:
        ENVS + "/python.yaml"
    shell:
        """
        python {params.script} --run_option coordinate \
            --inputFile1 {input.tmap} --inputFile2 {input.gtf_strg} > {output.tab} 2> {log}
        python {params.script} --run_option annotation \
            --inputFile1 {input.tmap} --inputFile2 {input.gtf_strg} > {output.gtf} 2>> {log}
        """


rule merged_trans_bdtM:
    input:
        expand("work/discovery/pics/{s}_bdtm.gtf", s=sorted(DISC_SAMPLES)),
    output:
        allpics="work/pics_annotation/allpics_annotation_bdtm.gtf",
        merged="work/pics_annotation/allpics_annotation_Merged_bdtm.bed",
    conda:
        ENVS + "/coco.yaml"
    log:
        "logs/discovery/merged_trans_bdtM.log",
    shell:
        """
        cat {input} > {output.allpics}
        # sort + strand-aware merge, then +1 on start (GTF-like 1-based coords)
        bedtools sort -i {output.allpics} \
            | bedtools merge -s -i stdin \
            | awk 'BEGIN{{OFS="\\t"}} {{print $1, $2+1, $3}}' > {output.merged} 2> {log}
        """


rule annotFile1_bdtM:
    input:
        bed="work/pics_annotation/allpics_annotation_Merged_bdtm.bed",
        gtf="work/pics_annotation/allpics_annotation_bdtm.gtf",
    output:
        "work/pics_annotation/pics_annotation_file_project1_bdtm.gtf",
    params:
        script=SCRIPTS + "/buildAnnotFileBdtm.py",
    log:
        "logs/discovery/annotFile1_bdtM.log",
    conda:
        ENVS + "/coco.yaml"
    shell:
        "python {params.script} -b {input.bed} -a {input.gtf} -o {output} &> {log}"


# ----------------------------------------------------------------- strgm ----
rule strg_merge:
    """StringTie --merge of all discovery samples (explicit, sorted file list)."""
    input:
        gtfs=expand("work/discovery/stringtie/{s}_strg.gtf", s=sorted(DISC_SAMPLES)),
        ref=REF_GTF,
    output:
        lst="work/discovery/strgm/stringtie_gtffile_list.txt",
        gtf="work/discovery/strgm/" + MERGED + "_strg.gtf",
    threads: 32
    resources:
        mem_mb=60000,
        runtime=240,
    log:
        "logs/discovery/strg_merge.log",
    conda:
        ENVS + "/stringtie.yaml"
    shell:
        """
        printf '%s\\n' {input.gtfs} > {output.lst}
        stringtie --merge -i -p {threads} -G {input.ref} -o {output.gtf} {output.lst} &> {log}
        """


rule gffcomp_strgM:
    input:
        gtf_ref=REF_GTF,
        gtf_strg="work/discovery/strgm/" + MERGED + "_strg.gtf",
    output:
        annotated="work/discovery/gffcomp_strgm/" + MERGED + ".annotated.gtf",
        tmap="work/discovery/strgm/" + MERGED + "." + MERGED + "_strg.gtf.tmap",
    params:
        prefix="work/discovery/gffcomp_strgm/" + MERGED,
    log:
        "logs/discovery/gffcomp_strgm.log",
    conda:
        ENVS + "/gffcompare.yaml"
    shell:
        "gffcompare -R -r {input.gtf_ref} -o {params.prefix} {input.gtf_strg} &> {log}"


rule target_trans_strgM:
    input:
        tmap=rules.gffcomp_strgM.output.tmap,
        gtf_strg=rules.gffcomp_strgM.input.gtf_strg,
    output:
        gtf="work/discovery/pics/" + MERGED + "_strgm.gtf",
        tab="work/discovery/pics/" + MERGED + "_strgm.tsv",
    params:
        script=SCRIPTS + "/getAnnotFileProject1.py",
    log:
        "logs/discovery/target_trans_strgM.log",
    conda:
        ENVS + "/python.yaml"
    shell:
        """
        python {params.script} --run_option coordinate \
            --inputFile1 {input.tmap} --inputFile2 {input.gtf_strg} > {output.tab} 2> {log}
        python {params.script} --run_option annotation \
            --inputFile1 {input.tmap} --inputFile2 {input.gtf_strg} > {output.gtf} 2>> {log}
        """


rule annotFile1_strgM:
    input:
        gtf=rules.target_trans_strgM.output.gtf,
    output:
        "work/pics_annotation/pics_annotation_file_project1_strgm.gtf",
    params:
        script=SCRIPTS + "/buildAnnotFileStrgm.py",
    log:
        "logs/discovery/annotFile1_strgM.log",
    conda:
        ENVS + "/python.yaml"
    shell:
        "python {params.script} -a {input.gtf} -o {output} &> {log}"
