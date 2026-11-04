# SnoDiscOv : snoRNA Discovery and Omics-based Validation 

A Snakemake-based bioinformatics pipeline for integrative transcriptomic discovery and characterization of human small nucleolar RNAs (snoRNAs).

**Associated publication:**

***Expanding the Landscape of Human snoRNAs through Integrative Transcriptomic Analyses***


snoDiscov runs in two stages:

1. **Discovery.** Candidate transcripts are assembled from total RNA-seq of human
   tissues with StringTie, using two complementary strategies (per-sample
   assemblies merged with bedtools, and `stringtie --merge`). Novel intergenic or
   intronic transfrags of 50–200 nt are kept if CoCo detects them at ≥ 1 TPM.
   → `results/gene_of_interest/gene_of_interest.{gtf,bed}`
2. **snoRNP dependence.** The candidates are added to the reference annotation and
   quantified in RNA-seq after depletion of snoRNP core proteins: DKC1 (H/ACA),
   NOP58 and SNU13/NHP2L1 (C/D), analysed with DESeq2. Their enrichment in SNU13
   RiPEN-seq is also measured.
   → `results/rbp_dep/`

Quality-control reports are in `results/qc/` and run provenance in
`results/provenance.yaml`; intermediate files are kept in `work/`.

The workflow is written in Snakemake and comes with a command-line wrapper, `snodiscov`.
It is tested with Snakemake 7.26 (the version used for the article) and Snakemake 8/9.

## Installation

```bash
git clone git@github.com:abthiaw/snoDiscOv.git && cd snodiscov
python -m venv ~/venvs/snodiscov && source ~/venvs/snodiscov/bin/activate
pip install .            # Snakemake 7.26+ (keeps an existing install), pandas, pyyaml
pip install .[slurm]     # optional: Snakemake >= 8 + SLURM executor plugin
```
You also need `conda` or `mamba` in your PATH (`conda config --set channel_priority strict` is recommended).
Every tool runs in a pinned conda environment (see [Software](#software)).

## Quick start

```bash
# 1. Create a run from a folder of FASTQ files (one sub-folder per dataset)
snodiscov init /path/to/fastq -o /scratch/snodiscov_run --resources-dir ~/projects/def-xxx/snodiscov_resources

# 2. On a node WITH internet (login / interactive node): download the reference
#    files from Zenodo (MD5-checked), install CoCo at a fixed commit, and build
#    all conda environments
snodiscov setup -d /scratch/snodiscov_run --conda-prefix ~/projects/def-xxx/snodiscov_envs

# 3. Run. Compute nodes need no internet access.
snodiscov run -d /scratch/snodiscov_run --profile slurm --account def-xxx --conda-prefix ~/projects/def-xxx/snodiscov_envs
snodiscov run -d /scratch/snodiscov_run -n                    # dry-run
snodiscov run -d RUN --profile local --cores 32               # single machine
snodiscov run -d RUN --profile slurm --account def-xxx -- --until coco_merge    # extra snakemake arguments
```
`snodiscov run` will not start if anything from `setup` is missing, so it never tries to download from a compute node.
Pass `--allow-downloads` if the node you run from does have internet access.

The built-in `slurm` and `local` profiles include one config for Snakemake 7 (`config.yaml`) and one for Snakemake 8+ (`config.v8+.yaml`).
Snakemake picks the right one automatically. Give your SLURM allocation with `--account def-xxx`, or set `SNODISCOV_SLURM_ACCOUNT`; the profile files do not need to be edited.

### FASTQ folder and manifest

```
fastq/
├── tissues/        Brain_1_R1.fastq.gz Brain_1_R2.fastq.gz ...      role: discovery
├── dkc1_kd/        DKC1_G_s_n1_S2_R1_001.fastq.gz ...               role: dkc1_depletion
├── nop58_kd/       NOP58_1_n2_R1.fastq.gz ...                       role: nop58_depletion
├── nhp2l1_kd/      NHP2L1_G_s_1_R1.fastq.gz LF_72h_1_R1.fastq.gz ... role: snu13_depletion
└── snu13_ripenseq/ SRR12998931_1.fastq.gz ...                       role: snu13_ripenseq
```
Each sub-folder name must match a key under `datasets:` in `config.yaml`; `--map OLD=NEW` renames a folder.
Read mates are detected from `_R1/_R2`, `_1/_2` or `_R1_001`. A file with no mate is treated as single-end.
Sample names can only contain `[A-Za-z0-9_-]`; other characters are replaced by `_`.

The `condition` column is filled from the regexes under `datasets.<name>.conditions`. **Check it before running DESeq2.**
For the depletion datasets the reference level must be called `control` (`deseq2_analysis_smk.R` relevels on it).
The DESeq2 design table (`results/rbp_dep/<target>_deseq2_design.tsv`, columns `samples` and `features`) is built from this column.

`manifest.tsv` columns: `sample dataset condition layout fq1 fq2 md5_fq1 md5_fq2` (add `--checksums` to fill the MD5 columns).
You can also write the manifest yourself. The config and the manifest are validated against JSON schemas (`workflow/schemas/`).

## Configuration

`config.yaml` (template: `snodiscov/config/config.yaml`):

* `resources`: every external file as a `zenodo: {record, file}`, `url` or `path` entry. The adapter files are shipped with the package (`snodiscov/data/Adapters/`, `package:` entries) and need nothing. Zenodo files are checked against the MD5 published by the Zenodo API, and each download leaves a `*.provenance.json` (URL, MD5, date).
* `tools`: git URL and **exact commit** of CoCo and pairedBamToBed12.
* `coverage_tracks` (default `false`): CoCo bedgraphs, bigWig files and bedgraph filtering for the discovery samples. These are long jobs that the main outputs do not need.
* `datasets`: role, adapters, cutadapt / Trimmomatic / STAR / CoCo arguments for each dataset. The defaults are the published parameters. An empty `cutadapt` or `trimmomatic` value skips that step.

## Outputs

```
results/
├── gene_of_interest/   main output 1 (+ SHA256SUMS)
├── rbp_dep/            main output 2 (+ SHA256SUMS)
├── qc/                 one MultiQC report per dataset: <dataset>_multiqc_report.html
└── provenance.yaml     versions, input checksums, tool commits and full config of the run
work/                   intermediate files (trimmed reads, BAM, CoCo counts, STAR indexes,
                        annotations...): kept so an interrupted run resumes without
                        recomputing; can be deleted once the results are archived
```

The two main outputs each have their own folder, with a `SHA256SUMS` file:

**`results/gene_of_interest/`**: expressed candidate transcripts (rule `coco_parsing`)

| File | Content |
|---|---|
| `gene_of_interest.gtf` | GTF (gene / transcript / exon) of the candidates with TPM ≥ 1 (`pics_gtffile1_project1.gtf` in the in-house pipeline) |
| `gene_of_interest.bed` | chr, start, end, gene_id, length, strand (`pics_bedfile1_project1.bed` in the in-house pipeline) |

**`results/rbp_dep/`**: snoRNP protein dependence (DESeq2 depletion analyses and RiPEN-seq)

| File | Content |
|---|---|
| `{dkc1,nop58,nhp2l1}_depletion_deseq2.tsv` | DESeq2 results (`deseq2_analysis_smk.R`) for the novel genes: Gene_Id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj |
| `{dkc1,nop58,nhp2l1}_deseq2_design.tsv` | DESeq2 design used (`samples`, `features`; reference level `control`) |
| `{dkc1,nop58,nhp2l1}_depletion_deseq2_Parsing.tsv` | per candidate (Id = chr:start-end): depletion call (`haca` for DKC1, `cd` for NOP58 / NHP2L1, else `negative`), l2fc, p-value |
| `snu13_ripenseq.tsv` | TPM of the candidates expressed (> 1 TPM) in SNU13 RiPEN-seq |
| `snu13_ripenseq_Parsing.tsv` | the same candidates labelled `cd` |

Other outputs:

| Path | Content |
|---|---|
| `results/qc/<dataset>_multiqc_report.html` | FastQC (raw + clean), cutadapt, Trimmomatic and STAR QC |
| `results/provenance.yaml` | versions, input checksums, tool commits and full config of the run |
| `work/<dataset>/coco_merge/*_merged_coco_{count,tpm}.tsv` | gene × sample matrices (DESeq2 input) |
| `work/annotation/*.gtf` | reference + candidate annotations (`annotFile1`, `annotFile2`) |
| `work/<discovery>/bigwig/*.bw` | coverage tracks (only with `coverage_tracks: true`) |

Runs made before this layout (everything under `results/`) can be converted
without recomputation: `bash tools/migrate_results_layout.sh <run_dir>`.

## Workflow overview

```
FASTQ ─ cutadapt ─ Trimmomatic ─ STAR(index: ref) ─ StringTie ┬ per-sample gffcompare + bedtools merge (bdtm)
                                                              └ StringTie --merge + gffcompare      (strgm)
        → PIC annotation → CoCo (discovery) → expressed PICs → annotFile2
        ─→ results/gene_of_interest/ (gene_of_interest.{gtf,bed})
KD / RiPEN-seq FASTQ ─ … ─ STAR(index: ref + PICs) ─ CoCo ─ DESeq2 (KD) / TPM (RiPEN-seq) ─→ results/rbp_dep/
```

Selection criteria (as implemented in the scripts):
* candidates (`coco_parsing`): StringTie transfrags of class `u` or `i` (gffcompare), 50–200 nt, TPM ≥ 1 (CoCo) in at least one discovery sample;
* depletion (`dep_parsing`, `picsParsing.deseq2Parsing`): DESeq2 p-value ≤ 0.05 and log2FC ≤ −0.5. Hits are labelled `haca` for DKC1 and `cd` for NOP58 and NHP2L1;
* SNU13 RiPEN-seq (`snu13_parsing`): TPM > 1 in at least one sample.

## Software

| Tool | Version | Env file |
|---|---|---|
| FastQC | 0.11.9 | `fastqc.yaml` |
| cutadapt | 4.6 | `cutadapt.yaml` |
| Trimmomatic | 0.36 | `trimmomatic.yaml` |
| STAR | 2.7.9a | `star.yaml` |
| StringTie | 2.1.5 | `stringtie.yaml` |
| gffcompare | 0.11.2 | `gffcompare.yaml` |
| bedtools | 2.31.0 | `bedtools.yaml`, `coco.yaml` |
| CoCo | git commit in `config.yaml` (+ subread 2.0.3, samtools 1.17) | `coco.yaml` |
| pairedBamToBed12 | git commit in `config.yaml` | built by `setup` |
| DESeq2 | 1.38.0 (R 4.2) | `DESeq2.yaml` |
| MultiQC | 1.14 † | `multiqc.yaml` |
| Python / pandas | 3.11.2 / 2.0.3 † | `python.yaml` |
| bedGraphToBigWig | UCSC 445 † | `ucsc.yaml` |

† not pinned in the in-house environments: confirm with `snodiscov export-envs` (below).

**Exact environments.** `snodiscov export-envs DIR -o env_pins` exports the package lists of the environments built by Snakemake.
`DIR` is a run directory, a `--conda-prefix` directory, or the original in-house pipeline directory: it reads `DIR/.snakemake/conda`.
It writes `<env>.linux-64.pin.txt` (`conda list --explicit`) and `<env>.lock.yaml` files.
Copy the pin files next to the env files in `snodiscov/workflow/envs/` and Snakemake uses them automatically, giving identical environments.

## Data availability

| Dataset | Role | Source/Accession number |
|---|---|---|
| Human tissues / TGIRTseq | discovery | PRJNA1478292, GSE126797, GSE157846|
| DKC1 knock-down (HCT116) | dkc1_depletion | GSE126797 |
| NOP58 knock-down | nop58_depletion | PRJNA1478292 |
| SNU13 (NHP2L1) knock-down | snu13_depletion | PRJNA1478292 |
| SNU13 RiPEN-seq | snu13_ripenseq | SRR12998931–SRR12998934, SRR15609476–SRR15609479 |
| Genome GRCh38 (Ensembl 109) | reference | Zenodo 8007660 |
| Annotation hg38 Ensembl 109 + Scott lab snoRNAs/tRNAs (2024) | reference | Zenodo 10713009 |

## FAIR summary

* **Findable**: versioned releases archived on Zenodo with a DOI (`.zenodo.json`, `CITATION.cff`). Every external file has a persistent identifier (Zenodo record + file name).
* **Accessible**: open-source code. Data come from public repositories (SRA/GEO, Zenodo) through `snodiscov setup`, which checks MD5 sums.
* **Interoperable**: standard formats (FASTQ, BAM, GTF, BED, bigWig, TSV), a JSON-schema-validated config and manifest, and a standard Snakemake layout that runs on SLURM or a single machine.
* **Reusable**: pinned conda environments and tool commits, documented parameters, explicit selection criteria, a `results/provenance.yaml` written after every run, a license, and tests (unit, DAG dry-run, end-to-end).

## Reproducibility notes (changes from the in-house version)

* StringTie gene_id suffixes no longer use `$SLURM_JOB_ID`. They are now deterministic (`geneid_suffix: index|sample`). Gene IDs of novel transcripts may therefore differ in name from the in-house version, but the loci are the same.
* File lists are explicit and sorted. Before, they came from `ls`/glob, so their content depended on what happened to be in the folder.
* All files the workflow reads or writes are declared, including the gffcompare `.tmap` files that used to be read implicitly.
* Paths are no longer hard-coded (`/home/...`, `../data`). Downloads are checksum-verified. External tools are pinned at a git commit. bedGraphToBigWig now comes from Bioconda instead of rsync.
* `slurmSubmit.py` with `--immediate-submit` is replaced by standard per-job submission. With immediate submission, local rules would start before their inputs exist on the cluster. SLURM requests now come from the `threads` and `resources` declared in each rule.
* Conda channels: `anaconda`/`defaults` (and the misspelled `default`) were removed. Only conda-forge and bioconda are used, with `nodefaults`.
* Script fixes. In each case the output is unchanged wherever the old code was correct:
  * `buildAnnotFileBdtm.py`: the first GTF line was written even when its start is not the start of a merged locus. That gave a gene with empty coordinates, which makes `bedtools sort` 2.31 crash.
  * `cocoParsingCount.py`: TPM is now read from the `tpm` column of the CoCo header. Before, it was `split()[4]`, and `split()` drops the empty `gene_name` field of novel genes, so the wrong column could be read, or the script could fail.
  * `deseq2_analysis_smk.R`: count columns are matched to the design table by sample name. Before, the script relied on `order()` and `arrange()` sorting the same way, which depends on the locale.
  * `cocoCount_merged.py` reads its files in sorted order. Shebangs are portable and the missing `import sys` was added.
  * `picsParsing.py`, `buildAnnotFile*.py`: `sys.stdout = open('/dev/stdout','w')` truncated the redirected logs (rule logs, SLURM output). It is replaced by `sys.stdout = sys.__stdout__`.
* Pipeline scope: the PAR-CLIP / eCLIP intersections and the final H/ACA / C/D integration (`PicsAnalysis.py`) were removed. The FBL and NOP56 knock-downs are not analysed. The main outputs are written to `results/gene_of_interest/` and `results/rbp_dep/`. Coverage tracks (`coco cb`, bigWig, bedgraph filter) are optional (`coverage_tracks`).
* Other bugs fixed: the SNU13 MultiQC used the tissue sample list; the second CoCo merge call overwrote the log of the first; Trimmomatic and FastQC thread counts did not match the SLURM allocation.

## Troubleshooting

**`LockException: Directory cannot be locked`**: a previous run was interrupted
(kill, disconnection, power loss) and left its lock, or another Snakemake is still
running on the same run directory. First check that nothing is running
(`tmux ls`, `ps -u $USER -f | grep snakemake`, `squeue -u $USER`), then remove the
lock without starting any job:
```bash
snodiscov run -d RUN -- --unlock
```
Files left incomplete by the interrupted run are redone automatically
(`rerun-incomplete` is set in the profiles).

**A job failed**: read its log in `RUN/logs/<dataset>/<step>/<sample>.log` (and
`RUN/logs/slurm/` on SLURM), fix the cause, then rerun the same `snodiscov run`
command. Only missing or failed outputs are recomputed.

**`CreateCondaEnvironmentException` / `LibMambaUnsatisfiableError`**: a conda
environment could not be solved. The message names the env file
(`snodiscov/workflow/envs/<tool>.yaml`) and the conflicting packages. Environments
are built by `snodiscov setup`, on a node with internet access: rerun it after any
change to an env file.

**`resource '<name>' ... not in Zenodo record`**: the record or the file name in
`RUN/config.yaml` (`resources: <name>`) does not match Zenodo. The message lists
the files the record actually contains. Edit the run config, not the template
installed with the package: `snodiscov init` copies the template once.

**`snodiscov: error: resources/tools are missing`**: `snodiscov setup -d RUN` has
not completed. Run it on a node with internet access before `snodiscov run`.

## Testing

```bash
pip install .[test]
pytest tests/                 # manifest builder
bash tests/dryrun.sh          # full DAG dry-run on empty FASTQs
bash tests/e2e/run_e2e.sh     # full run on synthetic data (needs stringtie, gffcompare, bedtools)
```
The end-to-end test runs every rule on a synthetic genome. StringTie `--merge`, gffcompare, bedtools and all project scripts run for real. FastQC, cutadapt, Trimmomatic, STAR, CoCo and DESeq2 are replaced by stand-ins (`tests/e2e/fakebin`). The test checks the export checksums, and that the loci designed as DKC1-, NOP58- and NHP2L1-dependent are detected.

## Citation and license

See `CITATION.cff` and `LICENSE`. Data: Zenodo records listed in `config.yaml`.
