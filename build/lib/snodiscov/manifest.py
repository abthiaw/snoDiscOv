"""Build a sample manifest (TSV) from a directory of FASTQ files.

Layout expected by default::

    FASTQ_DIR/
        tissues/        Brain_1_R1.fastq.gz  Brain_1_R2.fastq.gz ...
        dkc1_kd/        DKC1_G_s_n1_S2_R1_001.fastq.gz ...
        nop58_kd/       ...
        snu13_ripenseq/   SRR12998931_1.fastq.gz SRR12998931_2.fastq.gz ...

The first-level sub-directory gives the ``dataset`` column (it must match a key
of ``datasets:`` in config.yaml; use ``dataset_map`` to rename). FASTQ files
placed directly in FASTQ_DIR are assigned to ``default_dataset``.

Mates are detected with MATE_REGEX (``_R1``/``_R2``, ``_1``/``_2``, optional
Illumina ``_001`` suffix). A file without a mate is treated as single-end.
"""
from __future__ import annotations

import csv
import hashlib
import re
import sys
from collections import defaultdict
from pathlib import Path

FASTQ_EXT = re.compile(r"\.(fastq|fq)\.gz$")
MATE_REGEX = r"^(?P<sample>.+?)[._-]R?(?P<mate>[12])(?:_001)?\.(?:fastq|fq)\.gz$"
COLUMNS = ["sample", "dataset", "condition", "layout", "fq1", "fq2", "md5_fq1", "md5_fq2"]


class ManifestError(RuntimeError):
    pass


def _md5(path: Path, chunk: int = 1 << 22) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def infer_condition(sample: str, patterns: dict[str, str] | None) -> str:
    """Return the first condition whose regex matches the sample name."""
    for condition, regex in (patterns or {}).items():
        if re.search(regex, sample):
            return condition
    return "NA"


def scan_fastq_dir(
    fastq_dir: str | Path,
    default_dataset: str = "tissues",
    dataset_map: dict[str, str] | None = None,
    mate_regex: str = MATE_REGEX,
) -> dict[tuple[str, str], dict[str, Path]]:
    """Return {(dataset, sample): {"1": path, "2": path}}."""
    root = Path(fastq_dir).expanduser().resolve()
    if not root.is_dir():
        raise ManifestError(f"FASTQ directory not found: {root}")
    dataset_map = dataset_map or {}
    rx = re.compile(mate_regex)

    groups: dict[tuple[str, str], dict[str, Path]] = defaultdict(dict)
    files = sorted(p for p in root.rglob("*") if p.is_file() and FASTQ_EXT.search(p.name))
    if not files:
        raise ManifestError(f"No *.fastq.gz / *.fq.gz file found under {root}")

    for f in files:
        rel = f.relative_to(root)
        sub = rel.parts[0] if len(rel.parts) > 1 else None
        dataset = dataset_map.get(sub, sub) if sub else default_dataset
        m = rx.match(f.name)
        if m:
            key, mate = (dataset, m.group("sample")), m.group("mate")
        else:  # no mate token -> single-end, sample = file stem
            key, mate = (dataset, FASTQ_EXT.sub("", f.name)), "1"
        if mate in groups[key]:
            raise ManifestError(
                f"Duplicate mate {mate} for sample '{key[1]}' in dataset '{key[0]}':\n"
                f"  {groups[key][mate]}\n  {f}"
            )
        groups[key][mate] = f.resolve()

    # An orphan R2 (or an orphan '_1' that is really part of the name) -> single-end
    fixed: dict[tuple[str, str], dict[str, Path]] = {}
    for (dataset, sample), mates in groups.items():
        if set(mates) == {"1", "2"}:
            fixed[(dataset, sample)] = mates
        else:
            path = next(iter(mates.values()))
            stem = FASTQ_EXT.sub("", path.name)
            # keep the stripped name only when the token was an explicit R1
            name = sample if re.search(r"[._-]R1(?:_001)?$", stem) else stem
            print(f"[manifest] '{path.name}' has no mate -> single-end sample '{name}'",
                  file=sys.stderr)
            if (dataset, name) in fixed:
                raise ManifestError(f"Sample name collision: {dataset}/{name}")
            fixed[(dataset, name)] = {"1": path}
    return fixed


def build_manifest(
    fastq_dir: str | Path,
    output: str | Path,
    datasets_cfg: dict | None = None,
    default_dataset: str = "tissues",
    dataset_map: dict[str, str] | None = None,
    mate_regex: str = MATE_REGEX,
    checksums: bool = False,
) -> list[dict[str, str]]:
    """Scan *fastq_dir* and write a manifest TSV to *output*."""
    datasets_cfg = datasets_cfg or {}
    groups = scan_fastq_dir(fastq_dir, default_dataset, dataset_map, mate_regex)

    unknown = sorted({d for d, _ in groups} - set(datasets_cfg)) if datasets_cfg else []
    if unknown:
        raise ManifestError(
            "Dataset(s) not declared in config 'datasets:': " + ", ".join(unknown)
            + "\nRename the sub-directories, use --map SUBDIR=DATASET, or add them to the config."
        )

    rows = []
    for (dataset, sample), mates in sorted(groups.items()):
        clean = re.sub(r"[^A-Za-z0-9_-]", "_", sample)
        if clean != sample:
            print(f"[manifest] sample '{sample}' renamed '{clean}' (only [A-Za-z0-9_-] allowed)",
                  file=sys.stderr)
            sample = clean
        patterns = (datasets_cfg.get(dataset) or {}).get("conditions")
        pe = "2" in mates
        row = {
            "sample": sample,
            "dataset": dataset,
            "condition": infer_condition(sample, patterns),
            "layout": "PE" if pe else "SE",
            "fq1": str(mates["1"]),
            "fq2": str(mates["2"]) if pe else "",
            "md5_fq1": _md5(mates["1"]) if checksums else "",
            "md5_fq2": _md5(mates["2"]) if (checksums and pe) else "",
        }
        rows.append(row)

    # sanity: role datasets used for DESeq2 need >= 2 conditions
    for dataset, cfg in datasets_cfg.items():
        if (cfg or {}).get("role", "").endswith("_depletion"):
            conds = {r["condition"] for r in rows if r["dataset"] == dataset}
            if conds and ("NA" in conds or len(conds) < 2):
                print(f"[manifest] WARNING: dataset '{dataset}' conditions = {sorted(conds)}; "
                      "edit the 'condition' column before running DESeq2.", file=sys.stderr)

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows
