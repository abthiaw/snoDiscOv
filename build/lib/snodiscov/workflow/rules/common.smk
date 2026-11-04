# =============================================================================
#  Shared helpers: config/manifest validation, paths, input functions
# =============================================================================
import re
import sys

import pandas as pd
from snakemake.utils import validate

# make the snodiscov package importable (fetch helpers) when not pip-installed
sys.path.insert(0, os.path.abspath(os.path.join(workflow.basedir, "..", "..")))
from snodiscov.fetch import local_name  # noqa: E402

validate(config, schema="../schemas/config.schema.yaml")

SCRIPTS = os.path.join(workflow.basedir, "scripts")
PKG_DIR = os.path.abspath(os.path.join(workflow.basedir, ".."))  # snodiscov package
ENVS = "../envs"

# ---------------------------------------------------------------- manifest --
MANIFEST = pd.read_csv(config["manifest"], sep="\t", dtype=str, comment="#").fillna("")
validate(MANIFEST, schema="../schemas/manifest.schema.yaml")
MANIFEST = MANIFEST.set_index(["dataset", "sample"], drop=False).sort_index()

if MANIFEST.index.duplicated().any():
    dup = MANIFEST.index[MANIFEST.index.duplicated()].tolist()
    raise ValueError("Duplicated (dataset, sample) in manifest: " + str(dup))

_undeclared = sorted(set(MANIFEST["dataset"]) - set(config["datasets"]))
if _undeclared:
    raise ValueError("Manifest datasets missing from config 'datasets:': " + str(_undeclared))

DATASETS = sorted(set(MANIFEST["dataset"]))
DSCFG = {d: config["datasets"][d] for d in DATASETS}


def samples(dataset):
    return MANIFEST.loc[dataset, "sample"].tolist()


def layout(dataset, sample):
    return MANIFEST.loc[(dataset, sample), "layout"]


def is_pe(dataset, sample):
    return layout(dataset, sample) == "PE"


def mates(dataset, sample):
    return ["R1", "R2"] if is_pe(dataset, sample) else ["R1"]


def datasets_with_role(role):
    return [d for d in DATASETS if DSCFG[d]["role"] == role]


def dataset_for(role, required=True):
    found = datasets_with_role(role)
    if len(found) > 1:
        raise ValueError("More than one dataset with role '" + str(role) + "': " + str(found))
    if not found:
        if required:
            raise ValueError("No dataset with role '" + str(role) + "' in the manifest")
        return None
    return found[0]


DISCOVERY = dataset_for("discovery")
# Depletion / RiPEN-seq datasets are optional: each one present in the
# manifest adds its own outputs.
DKC1 = dataset_for("dkc1_depletion", required=False)
NOP58 = dataset_for("nop58_depletion", required=False)
SNU13 = dataset_for("snu13_ripenseq", required=False)
NHP2L1 = dataset_for("snu13_depletion", required=False)  # SNU13 (NHP2L1) knock-down
COVERAGE = bool(config.get("coverage_tracks", False))


def _alt(names):
    names = sorted(set(names), key=len, reverse=True)
    return "|".join(re.escape(n) for n in names) if names else "(?!x)x"


PE_SAMPLES = _alt(MANIFEST.loc[MANIFEST.layout == "PE", "sample"])
SE_SAMPLES = _alt(MANIFEST.loc[MANIFEST.layout == "SE", "sample"])


wildcard_constraints:
    dataset=_alt(DATASETS),
    sample=_alt(MANIFEST["sample"]),
    mate="R1|R2",
    index="discovery|validation",


# --------------------------------------------------------------- resources --
RESDIR = os.path.abspath(config["resources_dir"])
RES = {}  # resource key -> local path
for _key, _spec in config["resources"].items():
    if "package" in _spec:  # file shipped with snodiscov (e.g. data/Adapters/)
        RES[_key] = os.path.join(PKG_DIR, _spec["package"])
        if not os.path.isfile(RES[_key]):
            raise FileNotFoundError("resource '" + str(_key) + "': " + str(RES[_key]) + " not found in the package")
    elif "path" in _spec:  # local file provided by the user: used in place
        RES[_key] = os.path.abspath(os.path.expanduser(_spec["path"]))
    else:
        RES[_key] = os.path.join(RESDIR, "data", local_name(_key, _spec))

# adapter files must be FASTA ('>' headers): cutadapt rejects anything else
for _key in {c["adapters"] for c in DSCFG.values()}:
    if _key not in RES:
        raise ValueError("adapters '" + str(_key) + "' is not declared under resources:")
    if os.path.isfile(RES[_key]):
        with open(RES[_key]) as _fh:
            _first = next((l for l in _fh if l.strip()), "")
        if not _first.startswith(">"):
            raise ValueError("adapter file " + RES[_key] + " is not FASTA (first line: "
                             + repr(_first.strip()) + "); each sequence needs a '>name' header line")

TOOLS = os.path.join(RESDIR, "tools")
COCO_BIN = os.path.join(TOOLS, "coco", "bin")
COCO_DONE = os.path.join(TOOLS, "coco.done")
PB_BIN = os.path.join(TOOLS, "pairedBamToBed12", "bin")
PB_DONE = os.path.join(TOOLS, "pairedBamToBed12.done")

REF_FASTA = RES["genome_fasta"]
REF_GTF = RES["annotation_gtf"]


def adapters(wc):
    return RES[DSCFG[wc.dataset]["adapters"]]


def _fill(template, wc):
    return (template or "").replace("{adapters}", adapters(wc))


# ------------------------------------------------------- read-level helpers --
def raw_fastq(dataset, sample, mate):
    return "work/" + str(dataset) + "/00_raw/" + str(sample) + "_" + str(mate) + ".fastq.gz"


def clean_fastq(dataset, sample):
    """Final reads entering STAR, following the steps enabled for the dataset."""
    cfg = DSCFG[dataset]
    if cfg.get("trimmomatic"):
        step = "02_trimmomatic"
    elif cfg.get("cutadapt"):
        step = "01_cutadapt"
    else:
        step = "00_raw"
    return ["work/" + str(dataset) + "/" + str(step) + "/" + str(sample) + "_" + str(m) + ".fastq.gz" for m in mates(dataset, sample)]


def trimmomatic_input(wc):
    step = "01_cutadapt" if DSCFG[wc.dataset].get("cutadapt") else "00_raw"
    return ["work/" + str(wc.dataset) + "/" + str(step) + "/" + str(wc.sample) + "_" + str(m) + ".fastq.gz"
            for m in mates(wc.dataset, wc.sample)]


def index_for(dataset):
    """Discovery reads -> reference index; validation datasets -> index on the
    reference + discovered candidates (annotFile2), as in the original design."""
    return "discovery" if DSCFG[dataset]["role"] == "discovery" else "validation"


def star_index_marker(index):
    return "work/star_index/" + str(index) + "/chrNameLength.txt"


def bam(dataset, sample):
    return "work/" + str(dataset) + "/03_star/" + str(sample) + "/Aligned.sortedByCoord.out.bam"


def coco_gtf(dataset):
    if index_for(dataset) == "discovery":
        return "work/annotation/hg38_Ensembl_V109_Scottlab_2024_snoRNAs_tRNA_pics.correct_annotation.gtf"
    return "work/annotation/hg38_Ensembl_V109_Scottlab_2024_snoRNAs_pics_dep.correct_annotation.gtf"


def chrom_sizes():
    if "chrom_sizes" in RES:
        return RES["chrom_sizes"]
    return "work/star_index/discovery/chrNameLength_chr.txt"


# ------------------------------------------------- rule input functions --
D = DISCOVERY
DISC_SAMPLES = samples(D) if D else []
MERGED = config.get("merged_name", "human_tissus_samples")


def _geneid_suffix(wc):
    if config.get("geneid_suffix", "index") == "sample":
        return wc.sample
    return str(DISC_SAMPLES.index(wc.sample) + 1)


def _index_gtf(wc):
    if wc.index == "discovery":
        return REF_GTF
    return "work/annotation/hg38_Ensembl_V109_Scottlab_2024_snoRNAs_pics_dep.gtf"


def _fastqc_input(wc):
    if wc.stage == "raw":
        return raw_fastq(wc.dataset, wc.sample, wc.mate)
    clean = clean_fastq(wc.dataset, wc.sample)
    return clean[0] if wc.mate == "R1" else clean[1]


def _multiqc_input(wc):
    d = wc.dataset
    files = []
    for s in samples(d):
        for m in mates(d, s):
            files += expand("work/{d}/qc/fastqc_{st}/{s}_{m}_fastqc.zip",
                            d=d, st=["raw", "clean"], s=s, m=m)
        files.append("work/" + str(d) + "/03_star/" + str(s) + "/Log.final.out")
        if DSCFG[d].get("cutadapt"):
            files.append("logs/" + str(d) + "/cutadapt/" + str(s) + ".log")
        if DSCFG[d].get("trimmomatic"):
            files.append("logs/" + str(d) + "/trimmomatic/" + str(s) + ".log")
    return files


# depletion target (file prefix kept from the in-house pipeline) -> dataset role
DEP_TARGETS = {"dkc1": "dkc1_depletion", "nop58": "nop58_depletion", "nhp2l1": "snu13_depletion"}


def dep_dataset(target):
    return dataset_for(DEP_TARGETS[target])


# ------------------------------------------------------------------ targets --
def main_outputs():
    """The two main output directories {dir: [files]}.
    gene_of_interest: expressed candidates (rule coco_parsing);
    rbp_dep: DESeq2 depletion tables and SNU13 RiPEN-seq tables."""
    goi = ["results/gene_of_interest/gene_of_interest.gtf", "results/gene_of_interest/gene_of_interest.bed"]
    rbp = []
    for tag, d in (("dkc1", DKC1), ("nop58", NOP58), ("nhp2l1", NHP2L1)):
        if d:
            rbp += ["results/rbp_dep/" + str(tag) + "_" + str(x) for x in
                    ("deseq2_design.tsv", "depletion_deseq2.tsv", "depletion_deseq2_Parsing.tsv")]
    if SNU13:
        rbp += ["results/rbp_dep/snu13_ripenseq.tsv", "results/rbp_dep/snu13_ripenseq_Parsing.tsv"]
    return {"gene_of_interest": goi, "rbp_dep": rbp}


def all_targets():
    t = ["results/gene_of_interest/SHA256SUMS"]
    if any((DKC1, NOP58, NHP2L1, SNU13)):
        t.append("results/rbp_dep/SHA256SUMS")
    for d in DATASETS:
        t.append("results/qc/" + str(d) + "_multiqc_report.html")
        t += expand("work/{d}/coco_merge/{d}_merged_coco_{f}.tsv", d=d, f=["count", "tpm"])
    if COVERAGE:
        t += expand("work/{d}/bigwig/{s}_bigwig.bw", d=D, s=DISC_SAMPLES)
        t.append("work/discovery/bedgraph_filter/pics_bedgraph_filter_project1.tsv")
    return t


# --------------------------------------------------------------- provenance --
def _sha256(path, chunk=1 << 22):
    import hashlib

    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def write_provenance(out):
    import datetime
    import glob
    import json
    import platform

    import snakemake as _smk
    import yaml

    try:
        from snodiscov import __version__ as sv
    except Exception:  # noqa: BLE001
        sv = "unknown"
    tools = {}
    for t in ("coco", "pairedBamToBed12"):
        f = os.path.join(TOOLS, t, "COMMIT")
        tools[t] = {"repo": config["tools"][t]["repo"],
                    "commit": open(f).read().strip() if os.path.exists(f) else None}
    resources = {}
    for k, path in RES.items():
        prov = path + ".provenance.json"
        resources[k] = json.load(open(prov)) if os.path.exists(prov) else \
            {"path": path, "sha256": _sha256(path) if os.path.isfile(path) else None}
    envs = {os.path.basename(f): _sha256(f)
            for f in sorted(glob.glob(os.path.join(workflow.basedir, "envs", "*")))}
    info = {
        "snodiscov_version": sv,
        "snakemake_version": _smk.__version__,
        "python": platform.python_version(),
        "host": platform.node(),
        "finished": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "manifest": {"path": os.path.abspath(config["manifest"]),
                     "sha256": _sha256(config["manifest"])},
        "datasets": {d: len(samples(d)) for d in DATASETS},
        "tools": tools,
        "resources": resources,
        "conda_env_files_sha256": envs,
        "config": config,
    }
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        yaml.safe_dump(json.loads(json.dumps(info, default=str)), fh, sort_keys=False)
