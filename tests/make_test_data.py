"""Create an empty FASTQ tree mirroring the published datasets (for dry-runs/CI)."""
import gzip
import sys
from pathlib import Path

root = Path(sys.argv[1])
layout = {
    "tissues": (["Brain_1", "Brain_2", "Liver_1", "Liver_2", "HCT116_1"], "{s}_R{m}.fastq.gz"),
    "dkc1_kd": (["DKC1_G_s_n1_S2", "DKC1_G_s_n2_S5", "Lipofectamine_n1_S1", "Lipofectamine_n2_S4"],
                "{s}_R{m}_001.fastq.gz"),
    "nop58_kd": (["Allstar_n2", "Allstar_n3", "NOP58_1_n2", "NOP58_1_n3"], "{s}_R{m}.fastq.gz"),
    "snu13_ripenseq": (["SRR12998931", "SRR12998932"], "{s}_{m}.fastq.gz"),
    "nhp2l1_kd": (["LF_72h_1", "LF_72h_2", "NHP2L1_G_1_s_1", "NHP2L1_G_s_1"], "{s}_R{m}.fastq.gz"),
}
for ds, (samples, pat) in layout.items():
    d = root / ds
    d.mkdir(parents=True, exist_ok=True)
    for s in samples:
        for m in ([None] if "{m}" not in pat else [1, 2]):
            f = d / pat.format(s=s, m=m)
            with gzip.open(f, "wt") as fh:
                fh.write("@r\nACGT\n+\nIIII\n")
