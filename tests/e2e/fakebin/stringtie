#!/usr/bin/env python3
"""Stand-ins for heavy tools (end-to-end wiring test only)."""
import hashlib, os, re, shutil, subprocess, sys
from pathlib import Path

tool = Path(sys.argv[0]).name; a = sys.argv[1:]
LOCI = os.environ["E2E_LOCI"]

def opt(flag, n=1):
    i = a.index(flag); return a[i + 1] if n == 1 else a[i + 1:i + 1 + n]

if tool == "fastqc":
    out = opt("-o"); b = re.sub(r"\.(fastq|fq)\.gz$", "", Path(a[0]).name)
    for e in ("html", "zip"): Path(out, f"{b}_fastqc.{e}").write_text("fake")
elif tool == "cutadapt":
    outs = {opt("-o")} | ({opt("-p")} if "-p" in a else set())
    ins = [x for x in a if x.endswith(".fastq.gz") and x not in outs]
    shutil.copy(ins[0], opt("-o"))
    if "-p" in a: shutil.copy(ins[1], opt("-p"))
    print("This is cutadapt (fake)")
elif tool == "trimmomatic":
    files = [x for x in a if x.endswith(".fastq.gz")]
    if a[0] == "PE":
        for src, dst in [(0, 2), (0, 3), (1, 4), (1, 5)]: shutil.copy(files[src], files[dst])
    else: shutil.copy(files[0], files[1])
elif tool == "multiqc":
    Path(opt("-o")).mkdir(parents=True, exist_ok=True); Path(opt("-o"), opt("-n")).write_text("fake")
elif tool == "STAR":
    if opt("--runMode") == "genomeGenerate":
        Path(opt("--genomeDir"), "chrNameLength.txt").write_text("1\t200000\n")
    else:
        p = opt("--outFileNamePrefix")
        Path(p + "Aligned.sortedByCoord.out.bam").write_text(p)
        Path(p + "Log.final.out").write_text("fake")
elif tool == "stringtie":
    if "--merge" in a:
        sys.exit(subprocess.call([os.environ["REAL_STRINGTIE"]] + a))
    out, bam = opt("-o"), a[-1]; sample = Path(bam).parent.name
    h = int(hashlib.md5(sample.encode()).hexdigest(), 16)
    lines = ["# stringtie fake\n"]
    for k, li in enumerate(open(LOCI)):
        if (h >> k) % 4 == 0: continue                    # each sample misses some loci
        c, s, e, st = li.split(); j = (h >> (k + 8)) % 3  # small coordinate jitter
        s, e = int(s) + j, int(e) + j
        g = f'gene_id "MSTRG.{k+1}"; transcript_id "MSTRG.{k+1}.1";'
        lines.append(f'{c}\tStringTie\ttranscript\t{s}\t{e}\t1000\t{st}\t.\t{g} cov "50.0"; FPKM "100.0"; TPM "200.0";\n')
        lines.append(f'{c}\tStringTie\texon\t{s}\t{e}\t1000\t{st}\t.\t{g} exon_number "1"; cov "50.0";\n')
    Path(out).write_text("".join(lines))
elif tool == "coco":
    if a[0] == "ca":
        g = a[1]; shutil.copy(g, re.sub(r"\.gtf$", ".correct_annotation.gtf", g))
    elif a[0] == "cc":
        gtf, bam, out = a[-3:]; tag = Path(out).stem
        genes, start = [], {}
        for li in open(gtf):
            m = re.search(r'gene_id "([^"]+)"', li)
            if m and m.group(1) not in genes:
                genes.append(m.group(1)); start[m.group(1)] = int(li.split("\t")[3])
        loci = [int(l.split()[1]) for l in open(LOCI)]
        # DKC1 KD depletes loci 1-3 (H/ACA-like), NOP58 KD depletes loci 5-7 (C/D-like)
        # NHP2L1 (SNU13) KD depletes loci 6-8
        target = (loci[1:4] if tag.startswith("DKC1") else loci[5:8] if tag.startswith("NOP58")
                  else loci[6:9] if tag.startswith("NHP2L1") else [])
        dep = lambda g: any(abs(start[g] - t) <= 5 for t in target)
        # real CoCo layout: gene_id gene_name count cpm tpm (gene_name empty for novel genes)
        rows = ["gene_id\tgene_name\tcount\tcpm\ttpm\n"]
        for i, g in enumerate(genes):
            cnt = 400 // (8 if dep(g) else 1) + i
            name = "" if "STRG" in g else g + "_name"
            rows.append(f"{g}\t{name}\t{cnt}\t{cnt/2}\t{cnt/3:.2f}\n")
        Path(out).write_text("".join(rows))
    elif a[0] == "cb":
        out = a[-2]
        rows = ["track type=bedGraph\n"] + [f"chr1\t{i}\t{i+10}\t{5 + (i // 10) % 7}\n" for i in range(0, 60000, 10)]
        Path(out).write_text("".join(rows))
elif tool == "bedGraphToBigWig":
    Path(a[2]).write_text("fake")
elif tool == "Rscript":   # DESeq2 stand-in: parse the Snakemake R preamble
    src = Path(a[-1]).read_text()
    get = lambda k: re.search(rf'"{k}" = [\'"]([^\'"]+)', src).group(1)
    import csv
    cnt = list(csv.DictReader(open(get("Count")), delimiter="\t"))
    design = {r["samples"]: r["features"] for r in csv.DictReader(open(get("condition")), delimiter="\t")}
    with open(get("deseq2_out"), "w") as fh:
        fh.write("Gene_Id\tbaseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj\n")
        import math
        for r in cnt:
            if "STRG" not in r["gene_id"]: continue
            ctl = sum(float(r[s]) for s, f in design.items() if f == "control")
            kd = sum(float(r[s]) for s, f in design.items() if f != "control")
            l2 = math.log2((kd + 1) / (ctl + 1)); p = 0.001 if l2 < -0.5 else 0.6
            fh.write(f"{r['gene_id']}\t{ctl}\t{l2:.3f}\t0.1\t-5\t{p}\t{p}\n")
else:
    sys.exit(f"unknown fake tool {tool}")
