#!/usr/bin/env bash
# Full-DAG dry-run on empty FASTQs and placeholder references (used in CI).
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
T=$(mktemp -d)
python "$HERE/make_test_data.py" "$T/fastq"
mkdir -p "$T/refs" && for r in genome.fa ref.gtf a.fa; do touch "$T/refs/$r"; done
python - "$T" <<'PY'
import sys, yaml, pathlib
T = pathlib.Path(sys.argv[1]); R = T / "refs"
cfg = yaml.safe_load((pathlib.Path("snodiscov/config/config.yaml")).read_text())
cfg["resources"] = {k: (v if "package" in v else
                        {"path": str(R / ("genome.fa" if k == "genome_fasta" else "ref.gtf"))})
                    for k, v in cfg["resources"].items()}
for t in cfg["tools"].values(): t["ref"] = "HEAD"
(T / "base.yaml").write_text(yaml.safe_dump(cfg))
PY
snodiscov init "$T/fastq" -o "$T/run" -c "$T/base.yaml"
python -m snakemake --configfile "$T/run/config.yaml" --snakefile snodiscov/workflow/Snakefile \
    --directory "$T/run" -n --quiet rules --cores 1
