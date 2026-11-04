# =============================================================================
#  Read pre-processing, generic over datasets (PE or SE), driven by the manifest
#  00_raw (symlinks, normalised names) -> 01_cutadapt -> 02_trimmomatic -> QC
# =============================================================================


localrules:
    link_fastq,


rule link_fastq:
    """Expose manifest FASTQs under normalised names: {sample}_{R1,R2}.fastq.gz"""
    input:
        lambda wc: MANIFEST.loc[(wc.dataset, wc.sample), "fq1" if wc.mate == "R1" else "fq2"],
    output:
        "work/{dataset}/00_raw/{sample}_{mate}.fastq.gz",
    shell:
        "ln -sfn $(readlink -f {input}) {output}"


rule cutadapt_pe:
    input:
        r1="work/{dataset}/00_raw/{sample}_R1.fastq.gz",
        r2="work/{dataset}/00_raw/{sample}_R2.fastq.gz",
        adapters=adapters,
    output:
        r1="work/{dataset}/01_cutadapt/{sample}_R1.fastq.gz",
        r2="work/{dataset}/01_cutadapt/{sample}_R2.fastq.gz",
    wildcard_constraints:
        sample=PE_SAMPLES,
    params:
        args=lambda wc: _fill(DSCFG[wc.dataset]["cutadapt"], wc),
    threads: 8
    resources:
        mem_mb=20000,
        runtime=300,
    log:
        "logs/{dataset}/cutadapt/{sample}.log",
    conda:
        ENVS + "/cutadapt.yaml"
    shell:
        "cutadapt {params.args} -j {threads} -o {output.r1} -p {output.r2} "
        "{input.r1} {input.r2} &> {log}"


rule cutadapt_se:
    input:
        r1="work/{dataset}/00_raw/{sample}_R1.fastq.gz",
        adapters=adapters,
    output:
        r1="work/{dataset}/01_cutadapt/{sample}_R1.fastq.gz",
    wildcard_constraints:
        sample=SE_SAMPLES,
    params:
        args=lambda wc: _fill(DSCFG[wc.dataset]["cutadapt"], wc),
    threads: 8
    resources:
        mem_mb=20000,
        runtime=300,
    log:
        "logs/{dataset}/cutadapt/{sample}.log",
    conda:
        ENVS + "/cutadapt.yaml"
    shell:
        "cutadapt {params.args} -j {threads} -o {output.r1} {input.r1} &> {log}"


rule trimmomatic_pe:
    input:
        reads=trimmomatic_input,
        adapters=adapters,
    output:
        r1="work/{dataset}/02_trimmomatic/{sample}_R1.fastq.gz",
        r2="work/{dataset}/02_trimmomatic/{sample}_R2.fastq.gz",
        u1="work/{dataset}/02_trimmomatic/unpaired/{sample}_R1.fastq.gz",
        u2="work/{dataset}/02_trimmomatic/unpaired/{sample}_R2.fastq.gz",
    wildcard_constraints:
        sample=PE_SAMPLES,
    params:
        steps=lambda wc: _fill(DSCFG[wc.dataset]["trimmomatic"], wc),
    threads: 8
    resources:
        mem_mb=20000,
        runtime=120,
    log:
        "logs/{dataset}/trimmomatic/{sample}.log",
    conda:
        ENVS + "/trimmomatic.yaml"
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.reads} "
        "{output.r1} {output.u1} {output.r2} {output.u2} {params.steps} &> {log}"


rule trimmomatic_se:
    input:
        reads=trimmomatic_input,
        adapters=adapters,
    output:
        r1="work/{dataset}/02_trimmomatic/{sample}_R1.fastq.gz",
    wildcard_constraints:
        sample=SE_SAMPLES,
    params:
        steps=lambda wc: _fill(DSCFG[wc.dataset]["trimmomatic"], wc),
    threads: 8
    resources:
        mem_mb=20000,
        runtime=120,
    log:
        "logs/{dataset}/trimmomatic/{sample}.log",
    conda:
        ENVS + "/trimmomatic.yaml"
    shell:
        "trimmomatic SE -threads {threads} -phred33 {input.reads} {output.r1} "
        "{params.steps} &> {log}"


rule fastqc:
    input:
        _fastqc_input,
    output:
        html="work/{dataset}/qc/fastqc_{stage}/{sample}_{mate}_fastqc.html",
        zip="work/{dataset}/qc/fastqc_{stage}/{sample}_{mate}_fastqc.zip",
    wildcard_constraints:
        stage="raw|clean",
    params:
        outdir=lambda wc, output: os.path.dirname(output.html),
    threads: 4
    resources:
        mem_mb=8000,
        runtime=120,
    log:
        "logs/{dataset}/fastqc_{stage}/{sample}_{mate}.log",
    conda:
        ENVS + "/fastqc.yaml"
    shell:
        "fastqc {input} -f fastq -t {threads} -o {params.outdir} &> {log}"


rule multiqc:
    """One report per dataset: FastQC (raw + clean), cutadapt, Trimmomatic, STAR."""
    input:
        _multiqc_input,
    output:
        "results/qc/{dataset}_multiqc_report.html",
    params:
        outdir=lambda wc, output: os.path.dirname(output[0]),
    log:
        "logs/{dataset}/multiqc.log",
    conda:
        ENVS + "/multiqc.yaml"
    shell:
        "multiqc -f -o {params.outdir} -n {wildcards.dataset}_multiqc_report.html {input} &> {log}"
