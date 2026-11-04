#!/usr/bin/env bash
# End-to-end test: real StringTie --merge, gffcompare, bedtools and all
# project scripts; heavy steps (FastQC, cutadapt, Trimmomatic, STAR, CoCo,
# DESeq2) replaced by stand-ins in fakebin/. Needs stringtie, gffcompare and
# bedtools in PATH (REAL_STRINGTIE = path of the real stringtie binary).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); REPO=$(cd "$HERE/../.." && pwd)
T=${1:-$(mktemp -d)}
python "$HERE/make_data.py" "$T"
export E2E_LOCI="$T/loci.tsv" REAL_STRINGTIE=${REAL_STRINGTIE:-$(command -v stringtie)}
export PATH="$HERE/fakebin:$PATH"
R="$T/refs"
python - "$T" "$R" "$REPO" <<'PY'
import sys, yaml, pathlib
T, R, REPO = map(pathlib.Path, sys.argv[1:])
cfg = yaml.safe_load((REPO / "snodiscov/config/config.yaml").read_text())
m = {"genome_fasta": "genome.fa", "annotation_gtf": "ref.gtf"}
# adapters: the files shipped in snodiscov/data/Adapters (template `package:` entries)
cfg["resources"] = {k: v for k, v in cfg["resources"].items() if "package" in v}
cfg["resources"].update({k: {"path": str(R / v)} for k, v in m.items()})
for t in cfg["tools"].values(): t["ref"] = "HEAD"
(T / "base.yaml").write_text(yaml.safe_dump(cfg))
PY
snodiscov init "$T/fastq" -o "$T/run" -c "$T/base.yaml" --force
# fake tool installs (the real ones are git clones; not exercised here)
mkdir -p "$T/run/resources/tools/coco/bin" "$T/run/resources/tools/pairedBamToBed12/bin"
ln -sf "$HERE/fakebin/coco" "$T/run/resources/tools/coco/bin/coco"
touch "$T/run/resources/tools/coco.done" "$T/run/resources/tools/pairedBamToBed12.done"
python -m snakemake --configfile "$T/run/config.yaml" --snakefile "$REPO/snodiscov/workflow/Snakefile" \
    --directory "$T/run" --cores 4 ${SMK_EXTRA:-}
R="$T/run/results"
echo "== gene_of_interest"; ls "$R/gene_of_interest"
echo "== rbp_dep"; ls "$R/rbp_dep"
# expected from the synthetic design: DKC1 KD depletes loci 1-3, NOP58 KD 5-7, NHP2L1 KD 6-8
python - "$R" "$T/loci.tsv" <<'PY'
import sys, csv, hashlib
R, loci = sys.argv[1], [l.split() for l in open(sys.argv[2])]
for sub in ("gene_of_interest", "rbp_dep"):                 # checksums are right
    for line in open(f"{R}/{sub}/SHA256SUMS"):
        h, name = line.split()
        assert hashlib.sha256(open(f"{R}/{sub}/{name}", "rb").read()).hexdigest() == h, name
E = f"{R}/rbp_dep"
def hits(tag, box):
    rows = csv.DictReader(open(f"{E}/{tag}_depletion_deseq2_Parsing.tsv"), delimiter="\t")
    return sorted(int(r["Id"].split(":")[1].split("-")[0]) for r in rows if r["boxType"] == box)
near = lambda got, idx: all(any(abs(g - int(loci[i][1])) <= 5 for g in got) for i in idx)
dk, nop, nhp = hits("dkc1", "haca"), hits("nop58", "cd"), hits("nhp2l1", "cd")
assert near(dk, [1, 2, 3]) and len(dk) == 3, dk
assert near(nop, [5, 6, 7]) and len(nop) == 3, nop
assert near(nhp, [6, 7, 8]) and len(nhp) == 3, nhp
for t in ("dkc1", "nop58", "nhp2l1"):                       # raw DESeq2 table + design
    assert open(f"{E}/{t}_depletion_deseq2.tsv").readline().startswith("Gene_Id")
    assert open(f"{E}/{t}_deseq2_design.tsv").readline().split() == ["samples", "features"]
snu = list(csv.DictReader(open(f"{E}/snu13_ripenseq_Parsing.tsv"), delimiter="\t"))
assert snu and all(r["SNU13_ripenseq"] == "cd" for r in snu)
bed = [l.split("\t") for l in open(f"{R}/gene_of_interest/gene_of_interest.bed")]
assert len(bed) == len(loci) and all(len(b) == 6 for b in bed), len(bed)
print(f"E2E OK: {len(bed)} expressed candidates, DKC1-dependent {len(dk)}, NOP58-dependent "
      f"{len(nop)}, NHP2L1-dependent {len(nhp)}, SNU13 RiPEN {len(snu)}")
PY
