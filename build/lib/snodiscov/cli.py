"""snodiscov command line wrapper around the Snakemake workflow.

    snodiscov init     FASTQ_DIR -o RUN_DIR      manifest + config for a new run
    snodiscov manifest FASTQ_DIR -o manifest.tsv (re)build only the manifest
    snodiscov setup    -d RUN_DIR                downloads + tools + conda envs
                                                 (run on a node WITH internet)
    snodiscov run      -d RUN_DIR [--profile slurm] [-n] [-- <snakemake args>]
    snodiscov export-envs DIR -o env_pins       exact versions of built envs
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

from . import __version__
from .manifest import MATE_REGEX, ManifestError, build_manifest

PKG = Path(__file__).resolve().parent
SNAKEFILE = PKG / "workflow" / "Snakefile"
DEFAULT_CONFIG = PKG / "config" / "config.yaml"
PROFILES = PKG / "profiles"


def _smk_major() -> int:
    """Major version of the Snakemake installed next to this wrapper."""
    try:
        from importlib.metadata import version
        return int(version("snakemake").split(".")[0])
    except Exception:  # noqa: BLE001
        return 8


def _conda_flag() -> list[str]:
    # Snakemake 7: --use-conda ; Snakemake >= 8: --software-deployment-method conda
    return ["--use-conda"] if _smk_major() < 8 else ["--software-deployment-method", "conda"]


def _err(msg: str, code: int = 1):
    print(f"snodiscov: error: {msg}", file=sys.stderr)
    sys.exit(code)


def _parse_map(items: list[str] | None) -> dict[str, str]:
    out = {}
    for it in items or []:
        if "=" not in it:
            _err(f"--map expects SUBDIR=DATASET, got '{it}'")
        k, v = it.split("=", 1)
        out[k] = v
    return out


def _load_yaml(path: Path) -> dict:
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def _run_cfg(run_dir: Path) -> Path:
    cfg = run_dir / "config.yaml"
    if not cfg.exists():
        _err(f"{cfg} not found. Create the run with `snodiscov init` first.")
    return cfg


def _base_cmd(run_dir: Path, targets: list[str], conda_prefix: str | None) -> list[str]:
    # targets first: --configfile takes several values and would swallow them
    # (--configfile is also kept away from the end of the command so that
    # targets given after `--` by the user are not swallowed either)
    cmd = [sys.executable, "-m", "snakemake", *targets,
           "--configfile", str(_run_cfg(run_dir)),
           "--snakefile", str(SNAKEFILE),
           "--directory", str(run_dir)]
    if conda_prefix:
        cmd += ["--conda-prefix", conda_prefix]
    return cmd


def _snakemake(run_dir: Path, args: list[str], conda_prefix: str | None,
               targets: list[str] | None = None) -> int:
    cmd = _base_cmd(run_dir, targets or [], conda_prefix) + args
    print("+ " + " ".join(cmd), file=sys.stderr, flush=True)
    return subprocess.call(cmd)


def _profile_args(profile: str | None, account: str | None = None) -> list[str]:
    """--profile <dir>, plus the SLURM account injected into default-resources.

    Snakemake replaces (does not merge) the profile's default-resources when
    --default-resources is given, so the profile's defaults are re-read and
    passed together with slurm_account."""
    if not profile:
        return []
    p = Path(profile)
    if not p.exists():
        p = PROFILES / profile
    if not p.is_dir():
        _err(f"profile '{profile}' not found (built-in: {', '.join(x.name for x in PROFILES.iterdir())})")
    args = ["--profile", str(p)]
    cfg_file = p / "config.v8+.yaml" if _smk_major() >= 8 and (p / "config.v8+.yaml").exists() \
        else p / "config.yaml"
    cfg = _load_yaml(cfg_file) if cfg_file.exists() else {}
    uses_slurm = cfg.get("executor") == "slurm" or "slurm_account" in str(cfg.get("cluster", ""))
    if not uses_slurm:
        return args
    if not account:
        _err("a SLURM account is required with this profile: add --account def-xxx "
             "(or set SNODISCOV_SLURM_ACCOUNT)")
    res = cfg.get("default-resources") or []
    if isinstance(res, dict):
        res = [f"{k}={v}" for k, v in res.items()]
    res = [r for r in res if not str(r).startswith("slurm_account")]
    return args + ["--default-resources", f'slurm_account="{account}"', *map(str, res)]


# ------------------------------------------------------------------ commands --
def cmd_manifest(a) -> int:
    cfg = _load_yaml(Path(a.config)) if a.config else _load_yaml(DEFAULT_CONFIG)
    try:
        rows = build_manifest(a.fastq_dir, a.output, cfg.get("datasets", {}),
                              default_dataset=a.default_dataset,
                              dataset_map=_parse_map(a.map), mate_regex=a.pattern,
                              checksums=a.checksums)
    except ManifestError as e:
        _err(str(e))
    by_ds: dict[str, list] = {}
    for r in rows:
        by_ds.setdefault(r["dataset"], []).append(r)
    for ds, rs in sorted(by_ds.items()):
        lay = {r["layout"] for r in rs}
        conds = sorted({r["condition"] for r in rs})
        print(f"  {ds:<16} {len(rs):>3} samples  layout={'/'.join(sorted(lay))}  conditions={conds}")
    print(f"Manifest written: {a.output}  (review the 'condition' column before running)")
    return 0


def cmd_init(a) -> int:
    run_dir = Path(a.outdir).resolve()
    run_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = run_dir / "config.yaml"
    if cfg_path.exists() and not a.force:
        _err(f"{cfg_path} exists (use --force to overwrite)")
    src = Path(a.config) if a.config else DEFAULT_CONFIG
    text = src.read_text()
    cfg = yaml.safe_load(text)
    manifest = run_dir / "manifest.tsv"
    # keep comments: only rewrite the manifest / resources_dir lines
    lines = []
    for line in text.splitlines():
        if line.startswith("manifest:"):
            line = f"manifest: {manifest}"
        elif line.startswith("resources_dir:") and a.resources_dir:
            line = f"resources_dir: {Path(a.resources_dir).resolve()}"
        lines.append(line)
    cfg_path.write_text("\n".join(lines) + "\n")
    a.output, a.config = str(manifest), str(cfg_path)
    cmd_manifest(a)
    todo = [k for k, v in cfg.get("resources", {}).items()
            if "TODO" in str(v)] + [k for k, v in cfg.get("tools", {}).items() if "TODO" in str(v)]
    if todo:
        print(f"\nNOTE: fill the TODO entries in {cfg_path}: {', '.join(todo)}")
    print(f"\nNext:  snodiscov setup -d {run_dir}   (on a node with internet access)")
    return 0


def cmd_setup(a) -> int:
    run_dir = Path(a.directory).resolve()
    common = ["--cores", str(a.cores), *_conda_flag(), "--rerun-triggers", "mtime",
              "--rerun-incomplete", "--printshellcmds"]
    rc = _snakemake(run_dir, common, a.conda_prefix, targets=["setup"])
    if rc:
        return rc
    if not a.skip_envs:
        rc = _snakemake(run_dir, ["--conda-create-envs-only"] + common, a.conda_prefix)
    if rc == 0:
        print("\nSetup complete: references, CLIP data, tools and conda envs are local.\n"
              f"You can now run offline:  snodiscov run -d {run_dir} --profile slurm")
    return rc


def _setup_complete(run_dir: Path, conda_prefix: str | None) -> bool:
    cmd = _base_cmd(run_dir, ["setup"], conda_prefix) + [
        "--dry-run", "--rerun-triggers", "mtime", "--cores", "1"]
    out = subprocess.run(cmd, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stderr.write(out.stderr)
        _err("workflow could not be evaluated (see above)")
    txt = (out.stdout + out.stderr).lower()
    return "nothing to be done" in txt or not any(
        r in txt for r in ("fetch_resource", "install_coco", "install_pairedbamtobed12"))


def cmd_run(a) -> int:
    run_dir = Path(a.directory).resolve()
    if not a.allow_downloads and not _setup_complete(run_dir, a.conda_prefix):
        _err("resources/tools are missing. Run `snodiscov setup -d "
             f"{run_dir}` on a node with internet access first "
             "(or pass --allow-downloads if this node has internet).")
    args = _profile_args(a.profile, a.account)
    if a.cores:
        args += ["--cores", str(a.cores)]
    elif not a.profile:
        args += ["--cores", "all"]
    if a.dry_run:
        args += ["--dry-run"]
    args += _conda_flag()
    args += a.snakemake_args
    return _snakemake(run_dir, args, a.conda_prefix)


def _env_deps(path: Path) -> set[str]:
    try:
        deps = (_load_yaml(path) or {}).get("dependencies", []) or []
    except Exception:  # noqa: BLE001
        return set()
    return {str(d).split("=")[0].split()[0].lower() for d in deps if isinstance(d, str)}


def cmd_export_envs(a) -> int:
    """Export exact package lists of already-built Snakemake conda envs.

    Snakemake keeps, in .snakemake/conda/ (or --conda-prefix), a copy of every
    env file as <hash>.yaml next to the env prefix <hash>/. This works for the
    new runs AND for the in-house pipeline directory used for the article.
    """
    conda = shutil.which("mamba") or shutil.which("conda")
    if not conda:
        _err("conda/mamba not found in PATH")
    src = Path(a.conda_dir).expanduser().resolve()
    if (src / ".snakemake" / "conda").is_dir():
        src = src / ".snakemake" / "conda"
    out = Path(a.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    pkg_envs = {f.stem: _env_deps(f) for f in (PKG / "workflow" / "envs").glob("*.yaml")}
    platform = subprocess.run([conda, "info", "--json"], capture_output=True, text=True)
    try:
        import json
        platform = json.loads(platform.stdout).get("platform", "linux-64")
    except Exception:  # noqa: BLE001
        platform = "linux-64"
    n = 0
    for spec in sorted(src.glob("*.yaml")):
        prefix = spec.with_suffix("")
        if not prefix.is_dir():
            prefix = Path(str(spec)[:-5] + "_")      # Snakemake >= 8 suffixes '_'
            if not prefix.is_dir():
                continue
        deps = _env_deps(spec)
        name = next((k for k, v in pkg_envs.items() if v and v == deps), None) \
            or next((k for k, v in pkg_envs.items() if v and v <= deps), None) \
            or (sorted(deps)[0] if deps else spec.stem)
        with open(out / f"{name}.{platform}.pin.txt", "w") as fh:
            subprocess.run([conda, "list", "--explicit", "-p", str(prefix)], stdout=fh, check=True)
        with open(out / f"{name}.lock.yaml", "w") as fh:
            subprocess.run([conda, "env", "export", "--no-builds", "-p", str(prefix)], stdout=fh, check=True)
        print(f"  {spec.name:<40} -> {name}")
        n += 1
    if not n:
        _err(f"no built conda environment found in {src}")
    print(f"\n{n} environments exported to {out}.\nCopy the *.pin.txt files next to the "
          "env files in snodiscov/workflow/envs/ and update the versions in the .yaml files.")
    return 0


# --------------------------------------------------------------------- main --
def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="snodiscov", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"snodiscov {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    def manifest_opts(sp):
        sp.add_argument("fastq_dir", help="directory with FASTQ files (one sub-directory per dataset)")
        sp.add_argument("--map", action="append", metavar="SUBDIR=DATASET",
                        help="rename a sub-directory to a dataset name (repeatable)")
        sp.add_argument("--default-dataset", default="tissues",
                        help="dataset for FASTQs at the top level (default: tissues)")
        sp.add_argument("--pattern", default=MATE_REGEX,
                        help="regex with (?P<sample>) and (?P<mate>[12]) groups")
        sp.add_argument("--checksums", action="store_true", help="store MD5 of every FASTQ")

    sp = sub.add_parser("init", help="create a run directory (config + manifest)")
    manifest_opts(sp)
    sp.add_argument("-o", "--outdir", required=True, help="run directory")
    sp.add_argument("-c", "--config", help="base config (default: packaged template)")
    sp.add_argument("--resources-dir", help="shared directory for references/tools")
    sp.add_argument("--force", action="store_true")
    sp.set_defaults(func=cmd_init)

    sp = sub.add_parser("manifest", help="build only the manifest TSV")
    manifest_opts(sp)
    sp.add_argument("-o", "--output", required=True)
    sp.add_argument("-c", "--config", help="config providing datasets/conditions")
    sp.set_defaults(func=cmd_manifest)

    def run_opts(sp):
        sp.add_argument("-d", "--directory", required=True, help="run directory")
        sp.add_argument("--conda-prefix", default=os.environ.get("SNODISCOV_CONDA_PREFIX"),
                        help="shared conda env directory (env: SNODISCOV_CONDA_PREFIX)")

    sp = sub.add_parser("setup", help="download resources/tools and build conda envs (needs internet)")
    run_opts(sp)
    sp.add_argument("--cores", type=int, default=4)
    sp.add_argument("--skip-envs", action="store_true", help="do not pre-build conda envs")
    sp.set_defaults(func=cmd_setup)

    sp = sub.add_parser("run", help="run the workflow")
    run_opts(sp)
    sp.add_argument("--profile", help="'slurm', 'local' or a path to a Snakemake profile")
    sp.add_argument("--account", default=os.environ.get("SNODISCOV_SLURM_ACCOUNT"),
                    help="SLURM account (e.g. def-xxx), required with --profile slurm "
                         "(env: SNODISCOV_SLURM_ACCOUNT)")
    sp.add_argument("--cores", help="cores for local execution")
    sp.add_argument("-n", "--dry-run", action="store_true")
    sp.add_argument("--allow-downloads", action="store_true",
                    help="let local rules download missing resources during the run")
    sp.add_argument("snakemake_args", nargs=argparse.REMAINDER,
                    help="extra arguments passed to snakemake after '--'")
    sp.set_defaults(func=cmd_run)

    sp = sub.add_parser("export-envs", help="export exact versions of built conda envs (pin files)")
    sp.add_argument("conda_dir", help="a pipeline/run directory (uses .snakemake/conda) or a --conda-prefix dir")
    sp.add_argument("-o", "--output", default="env_pins")
    sp.set_defaults(func=cmd_export_envs)

    a = p.parse_args(argv)
    if getattr(a, "snakemake_args", None) and a.snakemake_args[0] == "--":
        a.snakemake_args = a.snakemake_args[1:]
    if shutil.which("conda") is None and shutil.which("mamba") is None and a.cmd in ("setup", "run"):
        print("snodiscov: warning: conda/mamba not found in PATH", file=sys.stderr)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
