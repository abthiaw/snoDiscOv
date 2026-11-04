# =============================================================================
#  Depletion analyses: DKC1, NOP58 and SNU13 (NHP2L1) knock-down (DESeq2),
#  and SNU13 RiPEN-seq (TPM), evaluated on the candidates of coco_parsing.
#  Main output 2: results/rbp_dep/
# =============================================================================
RBP = "results/rbp_dep"
DEP_GTF = "work/annotation/hg38_Ensembl_V109_Scottlab_2024_snoRNAs_pics_dep.correct_annotation.gtf"


wildcard_constraints:
    target="dkc1|nop58|nhp2l1",
    outdir="rbp_dep|gene_of_interest",


localrules:
    deseq2_conditions,
    checksums,


# ---------------------------------------------------------------- DESeq2 ----
rule deseq2_conditions:
    """Design table read by deseq2_analysis_smk.R: columns `samples` and
    `features` (the reference level must be named `control`). Built from the
    manifest, or copied from datasets.<name>.conditions_file if provided."""
    output:
        RBP + "/{target}_deseq2_design.tsv",
    params:
        dataset=lambda wc: dep_dataset(wc.target),
        user_file=lambda wc: DSCFG[dep_dataset(wc.target)].get("conditions_file", ""),
    log:
        "logs/rbp_dep/{target}_deseq2_design.log",
    run:
        import shutil

        if params.user_file:
            shutil.copyfile(params.user_file, output[0])
        else:
            tab = MANIFEST.loc[params.dataset, ["sample", "condition"]].rename(
                columns={"sample": "samples", "condition": "features"}
            )
            if "control" not in set(tab["features"]) or tab["features"].nunique() != 2:
                raise ValueError(
                    str(params.dataset) + ": manifest conditions must be 'control' and one "
                    "other level, got " + str(sorted(set(tab['features'])))
                )
            tab.sort_values("samples").to_csv(output[0], sep="\t", index=False)


rule deseq2:
    """deseq2_analysis_smk.R on the CoCo counts of a knock-down dataset."""
    input:
        Count=lambda wc: "work/{d}/coco_merge/{d}_merged_coco_count.tsv".format(d=dep_dataset(wc.target)),
        condition=RBP + "/{target}_deseq2_design.tsv",
    output:
        deseq2_out=RBP + "/{target}_depletion_deseq2.tsv",
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/rbp_dep/{target}_deseq2.log",
    conda:
        ENVS + "/DESeq2.yaml"
    script:
        "../scripts/deseq2_analysis_smk.R"


rule dep_parsing:
    """Candidate coordinates + depletion call (DepParsing.py): dkc1 hits are
    labelled 'haca', nop58 / nhp2l1 hits 'cd' (original depType labels)."""
    input:
        deseq2file=rules.deseq2.output.deseq2_out,
        snogtf=DEP_GTF,
    output:
        DepOut=RBP + "/{target}_depletion_deseq2_Parsing.tsv",
    params:
        depType="{target}_depletion_deseq2",
    conda:
        ENVS + "/python.yaml"
    log:
        "logs/rbp_dep/{target}_dep_parsing.log",
    script:
        "../scripts/DepParsing.py"


# ------------------------------------------------------------- RiPEN-seq ----
rule snu13_parsing:
    """Candidates with TPM > 1 in at least one SNU13 RiPEN-seq sample."""
    input:
        snu13_tpm=lambda wc: "work/{d}/coco_merge/{d}_merged_coco_tpm.tsv".format(
            d=dataset_for("snu13_ripenseq")
        ),
    output:
        snu13Out=RBP + "/snu13_ripenseq.tsv",
        snu13ParsOut=RBP + "/snu13_ripenseq_Parsing.tsv",
    params:
        script=SCRIPTS + "/snu13_snocand.py",
    conda:
        ENVS + "/python.yaml"
    log:
        "logs/rbp_dep/snu13_parsing.log",
    shell:
        """
        python {params.script} --input_count {input.snu13_tpm} \
            --output_count {output.snu13Out} --output_parsing {output.snu13ParsOut} &> {log}
        """


# ------------------------------------------------------------ checksums ----
rule checksums:
    """SHA256SUMS of each main output directory (integrity when copied/shared)."""
    input:
        lambda wc: main_outputs()[wc.outdir],
    output:
        "results/{outdir}/SHA256SUMS",
    log:
        "logs/{outdir}_checksums.log",
    run:
        with open(output[0], "w") as fh:
            for f in sorted(input):
                fh.write(str(_sha256(f)) + "  " + str(os.path.basename(f)) + "\n")
