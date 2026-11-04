"""Synthetic genome/annotation/StringTie-like data for the end-to-end test."""
import gzip, random, sys
from pathlib import Path

root = Path(sys.argv[1]); random.seed(1)
(root / "refs").mkdir(parents=True, exist_ok=True)
(root / "refs/genome.fa").write_text(">1\n" + "ACGT" * 50000 + "\n")
(root / "refs/ref.gtf").write_text(
    '1\tensembl\tgene\t1000\t5000\t.\t+\t.\tgene_id "ENSG01"; gene_biotype "protein_coding";\n'
    '1\tensembl\ttranscript\t1000\t5000\t.\t+\t.\tgene_id "ENSG01"; transcript_id "ENST01";\n'
    '1\tensembl\texon\t1000\t1500\t.\t+\t.\tgene_id "ENSG01"; transcript_id "ENST01"; exon_number "1";\n'
    '1\tensembl\texon\t4500\t5000\t.\t+\t.\tgene_id "ENSG01"; transcript_id "ENST01"; exon_number "2";\n')
# candidate loci: (start, length, strand) - one intronic (class i), rest intergenic (u)
LOCI = [(2500, 120, "+")] + [(10000 + 5000 * k, 70 + 10 * k, "+-"[k % 2]) for k in range(8)]
(root / "loci.tsv").write_text("".join(f"1\t{s}\t{s+l}\t{st}\n" for s, l, st in LOCI))

fq = root / "fastq"
def fastq(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh: fh.write("@r\nACGT\n+\nIIII\n")
for s in ["Brain_1", "Brain_2", "Liver_1"]:
    for m in (1, 2): fastq(fq / "tissues" / f"{s}_R{m}.fastq.gz")
for s in ["DKC1_G_s_n1_S2", "DKC1_G_s_n2_S5", "Lipofectamine_n1_S1", "Lipofectamine_n2_S4"]:
    for m in (1, 2): fastq(fq / "dkc1_kd" / f"{s}_R{m}_001.fastq.gz")
for s in ["Allstar_n2", "Allstar_n3", "NOP58_1_n2", "NOP58_1_n3"]:
    for m in (1, 2): fastq(fq / "nop58_kd" / f"{s}_R{m}.fastq.gz")
for s in ["LF_72h_1", "LF_72h_2", "NHP2L1_G_1_s_1", "NHP2L1_G_s_1"]:
    for m in (1, 2): fastq(fq / "nhp2l1_kd" / f"{s}_R{m}.fastq.gz")
for s in ["SRR12998931", "SRR12998932"]:
    for m in (1, 2): fastq(fq / "snu13_ripenseq" / f"{s}_{m}.fastq.gz")
