"""Download reference / CLIP resources (Zenodo or plain URL) with MD5 verification.

A resource spec (from config.yaml ``resources:``) is one of::

    {zenodo: {record: 8113986, file: DKC1_parclip_hg38.bed}, md5: <optional>, decompress: false}
    {url: https://..., md5: <optional>, decompress: true}
    {path: /local/file}                      # nothing to download

For Zenodo records the expected MD5 is taken from the Zenodo REST API when it
is not given in the config, so every file is checksum-verified.
Only the Python standard library is used: this runs on login / interactive
nodes, which are typically the only ones with internet access on HPC clusters.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.request
from pathlib import Path

ZENODO_API = os.environ.get("SNODISCOV_ZENODO_API", "https://zenodo.org/api")
UA = {"User-Agent": "snodiscov-fetch/1.0"}


def _log(msg: str) -> None:
    print(f"[fetch] {msg}", file=sys.stderr, flush=True)


def _get(url: str, timeout: int = 60):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout)


def zenodo_file_info(record: int | str, filename: str) -> tuple[str, str | None]:
    """Return (download_url, md5) for *filename* in Zenodo *record*."""
    with _get(f"{ZENODO_API}/records/{record}") as r:
        meta = json.load(r)
    for f in meta.get("files", []):
        if f.get("key") == filename:
            checksum = f.get("checksum", "")
            md5 = checksum.split(":", 1)[1] if checksum.startswith("md5:") else None
            url = f.get("links", {}).get("self") or \
                f"https://zenodo.org/records/{record}/files/{filename}?download=1"
            return url, md5
    available = ", ".join(f.get("key", "?") for f in meta.get("files", []))
    raise FileNotFoundError(f"'{filename}' not in Zenodo record {record} (has: {available})")


def md5sum(path: str | Path, chunk: int = 1 << 22) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _download(url: str, dest: Path, retries: int = 3) -> None:
    for attempt in range(1, retries + 1):
        try:
            _log(f"GET {url} -> {dest} (attempt {attempt})")
            with _get(url, timeout=300) as r, open(dest, "wb") as out:
                shutil.copyfileobj(r, out, length=1 << 22)
            return
        except Exception as exc:  # noqa: BLE001
            if attempt == retries:
                raise
            _log(f"failed ({exc}); retrying")
            time.sleep(5 * attempt)


def fetch_resource(spec: dict, output: str | Path) -> Path:
    """Materialise one resource at *output* (atomic, checksum-verified)."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    if "path" in spec:
        src = Path(spec["path"]).expanduser().resolve()
        if not src.exists():
            raise FileNotFoundError(src)
        if output.resolve() != src:
            output.unlink(missing_ok=True)
            output.symlink_to(src)
        return output

    md5 = spec.get("md5")
    if "zenodo" in spec:
        url, api_md5 = zenodo_file_info(spec["zenodo"]["record"], spec["zenodo"]["file"])
        md5 = md5 or api_md5
    elif "url" in spec:
        url = spec["url"]
    else:
        raise ValueError(f"Resource spec needs 'zenodo', 'url' or 'path': {spec}")

    raw = output.with_name(output.name + ".download")
    _download(url, raw)
    if md5:
        got = md5sum(raw)
        if got != md5:
            raw.unlink(missing_ok=True)
            raise ValueError(f"MD5 mismatch for {url}: expected {md5}, got {got}")
        _log(f"md5 OK ({md5})")
    else:
        _log("WARNING: no md5 available, file not verified")

    tmp = output.with_name(output.name + ".tmp")
    if spec.get("decompress"):
        with gzip.open(raw, "rb") as src, open(tmp, "wb") as dst:
            shutil.copyfileobj(src, dst, length=1 << 22)
        raw.unlink()
    else:
        raw.rename(tmp)
    tmp.rename(output)

    # provenance record next to the file
    prov = {"url": url, "md5_download": md5, "md5_final": md5sum(output),
            "retrieved": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    output.with_name(output.name + ".provenance.json").write_text(json.dumps(prov, indent=2))
    return output


def local_name(key: str, spec: dict) -> str:
    """File name used under resources_dir/data/ for a resource."""
    if "path" in spec or "package" in spec:
        name = Path(spec.get("path") or spec["package"]).name
    elif "zenodo" in spec:
        name = spec["zenodo"]["file"]
    else:
        name = spec["url"].split("?")[0].rstrip("/").split("/")[-1]
    if spec.get("rename"):
        name = spec["rename"]
    if spec.get("decompress") and name.endswith(".gz"):
        name = name[:-3]
    return name
