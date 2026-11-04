Adapter sequences used by cutadapt (`-a/-A file:`) and Trimmomatic (`ILLUMINACLIP`).

| File | Datasets |
|---|---|
| Adapters-PE_NextSeq.fa | tissues, NOP58 KD (cutadapt + Trimmomatic), DKC1 KD (Trimmomatic) |
| Adapters-PE_TruSeq.fa | SNU13 RiPEN-seq (RA5 / RA3 small-RNA adapters) |
| Adapters-PE_RBP_Toronto.fa | SNU13 (NHP2L1) KD |

Adapters-PE_RBP_Toronto.fa: the '>' of the two FASTA headers was restored
(the shared copy had `adapter_R1` / `adapter_R2` without '>', which cutadapt
rejects); sequences are unchanged.
