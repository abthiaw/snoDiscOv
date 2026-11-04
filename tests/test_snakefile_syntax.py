"""Snakemake 7.x re-tokenises Snakefiles; under Python >= 3.12 it corrupts
f-strings (inserts spaces, drops '{{' escapes). Keep rule files f-string free."""
import re
from pathlib import Path

WF = Path(__file__).resolve().parent.parent / "snodiscov" / "workflow"


def test_no_fstrings_in_rule_files():
    bad = []
    for f in [WF / "Snakefile", *sorted((WF / "rules").glob("*.smk"))]:
        for n, line in enumerate(f.read_text().splitlines(), 1):
            code = line.split("#", 1)[0]
            if re.search(r"(?<![A-Za-z0-9_])[fF][rR]?[\"']", code):
                bad.append(f"{f.name}:{n}: {line.strip()}")
    assert not bad, "f-strings break Snakemake 7 on Python 3.12:\n" + "\n".join(bad)
