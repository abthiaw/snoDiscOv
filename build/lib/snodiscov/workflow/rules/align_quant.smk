# =============================================================================
#  Alignment (STAR), quantification and coverage (CoCo), generic over datasets
#  Two STAR indexes:
#    discovery  : genome + reference GTF                  (tissues)
#    validation : genome + reference GTF + candidates     (KD / RiPEN-seq)
# =============================================================================


localrules:
    chrnamelength_chr,


rule star_index:
    input:
        fasta=REF_FASTA,
        gtf=_index_gtf,
    output:
        chrNameLength="work/star_index/{index}/chrNameLength.txt",
    params:
        outdir=lambda wc, output: os.path.dirname(output.chrNameLength),
        overhang=config["star_sjdb_overhang"],
    threads: 32
    resources:
        mem_mb=120000,
        runtime=240,
    log:
        "logs/star_index/{index}.log",
    conda:
        ENVS + "/star.yaml"
    shell:
        """
        rm -rf {params.outdir}/_STARtmp
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.outdir} \
            --outTmpDir {params.outdir}/_STARtmp \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            &> {log}
        """


rule chrnamelength_chr:
    """chrNameLength with a 'chr' prefix (chromosome sizes for bedGraphToBigWig)."""
    input:
        "work/star_index/{index}/chrNameLength.txt",
    output:
        "work/star_index/{index}/chrNameLength_chr.txt",
    conda:
        ENVS + "/coreutils.yaml"
    shell:
        """awk 'BEGIN{{OFS="\\t"}} {{print "chr"$1, $2}}' {input} > {output}"""


rule star_align:
    input:
        reads=lambda wc: clean_fastq(wc.dataset, wc.sample),
        index=lambda wc: star_index_marker(index_for(wc.dataset)),
    output:
        bam="work/{dataset}/03_star/{sample}/Aligned.sortedByCoord.out.bam",
        log_final="work/{dataset}/03_star/{sample}/Log.final.out",
    params:
        prefix=lambda wc, output: os.path.dirname(output.bam) + "/",
        index_dir=lambda wc, input: os.path.dirname(input.index),
        args=lambda wc: DSCFG[wc.dataset]["star"],
    threads: 32
    resources:
        mem_mb=100000,
        runtime=240,
    log:
        "logs/{dataset}/star/{sample}.log",
    conda:
        ENVS + "/star.yaml"
    shell:
        """
        STAR --runMode alignReads \
            --genomeDir {params.index_dir} \
            --readFilesIn {input.reads} \
            --runThreadN {threads} \
            --readFilesCommand zcat \
            --outStd Log \
            --outFileNamePrefix {params.prefix} \
            {params.args} \
            &> {log}
        """


rule coco_cc:
    """CoCo correct_count (strand-specific; --paired for PE libraries)."""
    input:
        gtf=lambda wc: coco_gtf(wc.dataset),
        bam=lambda wc: bam(wc.dataset, wc.sample),
        coco=COCO_DONE,
    output:
        tab="work/{dataset}/coco_cc/{sample}.tsv",
    params:
        coco_bin=COCO_BIN,
        args=lambda wc: DSCFG[wc.dataset].get("coco_cc", "--strand 1"),
        paired=lambda wc: "--paired" if is_pe(wc.dataset, wc.sample) else "",
    threads: 32
    resources:
        mem_mb=60000,
        runtime=1440,
    log:
        "logs/{dataset}/coco_cc/{sample}.log",
    conda:
        ENVS + "/coco.yaml"
    shell:
        """
        export PATH={params.coco_bin}:$PATH
        coco cc {params.args} --thread {threads} {params.paired} \
            {input.gtf} {input.bam} {output.tab} &> {log}
        """


rule coco_merge:
    """Merge per-sample CoCo tables into count and TPM matrices."""
    input:
        lambda wc: expand("work/{d}/coco_cc/{s}.tsv", d=wc.dataset, s=samples(wc.dataset)),
    output:
        counts="work/{dataset}/coco_merge/{dataset}_merged_coco_count.tsv",
        tpm="work/{dataset}/coco_merge/{dataset}_merged_coco_tpm.tsv",
    params:
        indir="work/{dataset}/coco_cc/",
        script=SCRIPTS + "/cocoCount_merged.py",
    log:
        "logs/{dataset}/coco_merge.log",
    conda:
        ENVS + "/python.yaml"
    shell:
        """
        (python {params.script} -i {params.indir} -o {output.counts} -f count
         python {params.script} -i {params.indir} -o {output.tpm} -f tpm) &> {log}
        """


rule coco_cb:
    """CoCo correct_bedgraph"""
    input:
        bam=lambda wc: bam(wc.dataset, wc.sample),
        chrNameLength=lambda wc: star_index_marker(index_for(wc.dataset)),
        coco=COCO_DONE,
        pb=PB_DONE,
    output:
        bg=temp("work/{dataset}/coco_cb/{sample}_unsorted.bedgraph"),
    params:
        coco_bin=COCO_BIN,
        pb_bin=PB_BIN,
    threads: 32
    resources:
        mem_mb=60000,
        runtime=720,
    log:
        "logs/{dataset}/coco_cb/{sample}.log",
    conda:
        ENVS + "/coco.yaml"
    shell:
        """
        export PATH={params.coco_bin}:{params.pb_bin}:$PATH
        coco cb -u -t {threads} -c 2500000 \
            {input.bam} {output.bg} {input.chrNameLength} &> {log}
        """


rule sort_bg:
    input:
        "work/{dataset}/coco_cb/{sample}_unsorted.bedgraph",
    output:
        "work/{dataset}/coco_cb/{sample}_sorted.bedgraph",
    resources:
        mem_mb=16000,
        runtime=120,
    conda:
        ENVS + "/coreutils.yaml"
    shell:
        """
        awk 'NR!=1' {input} \
            | LC_ALL=C sort -k1,1 -k2,2n -S 4G -T {resources.tmpdir} \
            | sed 's/chrM/chrMT/g' > {output}
        """


rule bedgraph_to_bigwig:
    input:
        bg="work/{dataset}/coco_cb/{sample}_sorted.bedgraph",
        sizes=chrom_sizes(),
    output:
        "work/{dataset}/bigwig/{sample}_bigwig.bw",
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/{dataset}/bigwig/{sample}.log",
    conda:
        ENVS + "/ucsc.yaml"
    shell:
        "bedGraphToBigWig {input.bg} {input.sizes} {output} &> {log}"
