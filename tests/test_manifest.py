import csv
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from snodiscov.manifest import ManifestError, build_manifest

HERE = Path(__file__).parent
CFG = yaml.safe_load((HERE.parent / "snodiscov" / "config" / "config.yaml").read_text())


@pytest.fixture
def fastq_dir(tmp_path):
    subprocess.check_call([sys.executable, str(HERE / "make_test_data.py"), str(tmp_path / "fq")])
    return tmp_path / "fq"


def read(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_layouts_and_names(fastq_dir, tmp_path):
    out = tmp_path / "m.tsv"
    build_manifest(fastq_dir, out, CFG["datasets"])
    rows = {(r["dataset"], r["sample"]): r for r in read(out)}
    assert rows[("dkc1_kd", "DKC1_G_s_n1_S2")]["layout"] == "PE"           # _R1_001
    assert rows[("snu13_ripenseq", "SRR12998931")]["fq2"].endswith("_2.fastq.gz")  # _1/_2
    assert rows[("tissues", "Brain_1")]["layout"] == "PE"
    assert rows[("dkc1_kd", "Lipofectamine_n1_S1")]["condition"] == "control"
    assert rows[("nop58_kd", "NOP58_1_n2")]["condition"] == "knockdown"
    assert rows[("nhp2l1_kd", "LF_72h_1")]["condition"] == "control"
    assert rows[("nhp2l1_kd", "NHP2L1_G_s_1")]["condition"] == "knockdown"


def test_unknown_dataset(fastq_dir, tmp_path):
    (fastq_dir / "weird").mkdir()
    (fastq_dir / "weird" / "x_R1.fastq.gz").write_bytes(b"")
    with pytest.raises(ManifestError):
        build_manifest(fastq_dir, tmp_path / "m.tsv", CFG["datasets"])


def test_dataset_map(fastq_dir, tmp_path):
    (fastq_dir / "tissues").rename(fastq_dir / "TOTAL_RNA")
    build_manifest(fastq_dir, tmp_path / "m.tsv", CFG["datasets"], dataset_map={"TOTAL_RNA": "tissues"})
    assert any(r["dataset"] == "tissues" for r in read(tmp_path / "m.tsv"))


def test_single_end_detection(tmp_path):
    """Single-end libraries (no mate file) stay supported for future datasets."""
    from snodiscov.manifest import scan_fastq_dir
    d = tmp_path / "fq" / "tissues"; d.mkdir(parents=True)
    for n in ["A_rep1.fastq.gz", "B_R1.fastq.gz", "C_R1.fastq.gz", "C_R2.fastq.gz"]:
        (d / n).write_bytes(b"")
    g = scan_fastq_dir(tmp_path / "fq")
    assert set(g[("tissues", "A_rep1")]) == {"1"}          # no mate token
    assert set(g[("tissues", "B")]) == {"1"}               # explicit R1, no R2
    assert set(g[("tissues", "C")]) == {"1", "2"}
