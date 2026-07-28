#!/usr/bin/env python3
"""Download all public SRA and GEO data linked to an NCBI BioProject."""
# Usage: DownloadBioproject.py PRJNA1147209 --workers 6 --output-root ./

from __future__ import annotations

import argparse
import csv
import hashlib
import html.parser
import json
import logging
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote, urljoin


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
SRA_LOCATOR = "https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"
BIOPROJECT_PATTERN = re.compile(r"^PRJNA[0-9]+$")
GEO_PATTERN = re.compile(r"^GSE[0-9]+$")


@dataclass(frozen=True)
class SraFile:
    """A downloadable full-quality SRA archive."""

    run_accession: str
    biosample: str
    sample_name: str
    size_bytes: int
    md5: str
    url: str


@dataclass(frozen=True)
class RemoteFile:
    """A remote file and its available validation metadata."""

    url: str
    output_path: Path
    size_bytes: int | None = None
    md5: str | None = None


class LinkParser(html.parser.HTMLParser):
    """Collect href attributes from an HTTP directory listing."""

    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(
        self, tag: str, attrs: list[tuple[str, str | None]]
    ) -> None:
        if tag != "a":
            return
        for name, value in attrs:
            if name == "href" and value:
                self.links.append(value)


def parse_args() -> argparse.Namespace:
    """Parse command-line parameters."""
    parser = argparse.ArgumentParser(
        description=(
            "Download full-quality SRA archives and all GEO files linked to "
            "an NCBI BioProject."
        )
    )
    parser.add_argument("bioproject", help="NCBI accession, for example PRJNA1147209")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("downloads"),
        help="Parent output directory (default: downloads)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Concurrent large-file downloads (default: 4)",
    )
    return parser.parse_args()


def run_curl(
    url: str,
    *,
    output_path: Path | None = None,
    resume: bool = False,
    capture: bool = False,
) -> str:
    """Run curl with retries and optional resumable output."""
    command = [
        "curl.exe" if sys.platform == "win32" else "curl",
        "-L",
        "--fail",
        "--show-error",
        "--retry",
        "20",
        "--retry-all-errors",
        "--retry-delay",
        "10",
        "--connect-timeout",
        "30",
    ]
    if capture:
        command.append("--silent")
    if resume:
        command.extend(["--continue-at", "-"])
    if output_path is not None:
        command.extend(["--output", str(output_path)])
    command.append(url)
    result = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE if capture else None,
    )
    if result.returncode != 0:
        raise RuntimeError(f"curl failed with exit code {result.returncode}: {url}")
    return result.stdout if capture else ""


def query_json(url: str) -> dict[str, Any]:
    """Retrieve and decode a JSON response."""
    value = json.loads(run_curl(url, capture=True))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object from {url}")
    return value


def write_text_new(path: Path, text: str) -> None:
    """Write a new metadata file without replacing an existing file."""
    if path.exists():
        if path.stat().st_size == 0:
            raise ValueError(f"Existing metadata file is empty: {path}")
        return

    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(text)


def md5sum(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    """Calculate an MD5 digest without loading the file into memory."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def verify_file(remote: RemoteFile) -> bool:
    """Return whether an existing final file matches all available metadata."""
    if not remote.output_path.is_file():
        return False
    if remote.size_bytes is not None:
        if remote.output_path.stat().st_size != remote.size_bytes:
            raise ValueError(f"Existing final file has a wrong size: {remote.output_path}")
    if remote.md5 is not None and md5sum(remote.output_path) != remote.md5.lower():
        raise ValueError(f"Existing final file failed MD5: {remote.output_path}")
    return True


def download_file(remote: RemoteFile) -> None:
    """Resume, validate, and atomically finalize one download."""
    if verify_file(remote):
        logging.info("SKIP verified file: %s", remote.output_path)
        return
    partial_path = remote.output_path.with_name(remote.output_path.name + ".part")
    logging.info("DOWNLOAD start: %s", remote.url)
    run_curl(remote.url, output_path=partial_path, resume=True)
    if remote.size_bytes is not None and partial_path.stat().st_size != remote.size_bytes:
        raise ValueError(f"Downloaded file has a wrong size: {partial_path}")
    if remote.md5 is not None and md5sum(partial_path) != remote.md5.lower():
        raise ValueError(f"Downloaded file failed MD5: {partial_path}")
    partial_path.replace(remote.output_path)
    logging.info("DOWNLOAD verified: %s", remote.output_path)


def batched(values: list[str], size: int) -> Iterable[list[str]]:
    """Yield fixed-size list batches."""
    for start in range(0, len(values), size):
        yield values[start : start + size]


def fetch_bioproject_xml(accession: str, metadata_dir: Path) -> tuple[str, list[str]]:
    """Download BioProject XML and return its numeric ID and linked GEO series."""
    search_url = (
        f"{EUTILS}/esearch.fcgi?db=bioproject"
        f"&term={quote(accession + '[Project Accession]')}&retmode=json&retmax=2"
    )
    ids = query_json(search_url)["esearchresult"]["idlist"]
    if len(ids) != 1:
        raise ValueError(f"Expected one BioProject record for {accession}, found {len(ids)}")
    project_id = str(ids[0])
    xml_url = f"{EUTILS}/efetch.fcgi?db=bioproject&id={project_id}&retmode=xml"
    xml_text = run_curl(xml_url, capture=True)
    write_text_new(metadata_dir / "bioproject.xml", xml_text)
    root = ET.fromstring(xml_text)
    geo_accessions = sorted(
        {
            (node.text or "").strip()
            for node in root.iter()
            if GEO_PATTERN.fullmatch((node.text or "").strip())
        }
    )
    return project_id, geo_accessions


def fetch_runinfo(accession: str, metadata_dir: Path) -> list[dict[str, str]]:
    """Download and validate SRA RunInfo rows associated with the project."""
    search_url = (
        f"{EUTILS}/esearch.fcgi?db=sra"
        f"&term={quote(accession + '[BioProject]')}&retmode=json&retmax=100000"
    )
    ids = [str(value) for value in query_json(search_url)["esearchresult"]["idlist"]]
    if not ids:
        raise ValueError(f"No public SRA runs found for {accession}")
    csv_parts: list[str] = []
    for id_batch in batched(ids, 200):
        fetch_url = (
            f"{EUTILS}/efetch.fcgi?db=sra&id={','.join(id_batch)}"
            "&rettype=runinfo&retmode=text"
        )
        csv_parts.append(run_curl(fetch_url, capture=True))
    rows: list[dict[str, str]] = []
    fieldnames: list[str] | None = None
    for csv_text in csv_parts:
        reader = csv.DictReader(csv_text.splitlines())
        if reader.fieldnames is None:
            raise ValueError("SRA RunInfo response has no header")
        fieldnames = reader.fieldnames
        rows.extend(dict(row) for row in reader)
    if len(rows) != len(ids):
        raise ValueError(f"Expected {len(ids)} RunInfo rows, received {len(rows)}")
    if {row["BioProject"] for row in rows} != {accession}:
        raise ValueError("SRA RunInfo contains unexpected BioProject accessions")
    runinfo_path = metadata_dir / "sra_runinfo.csv"
    if not runinfo_path.exists():
        with runinfo_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
    return rows


def locate_sra(row: dict[str, str]) -> SraFile:
    """Resolve one run to the full-quality normalized SRA object."""
    response = query_json(f"{SRA_LOCATOR}?acc={quote(row['Run'])}")
    results = response.get("result", [])
    if not results or results[0].get("status") != 200:
        raise ValueError(f"NCBI locator failed for {row['Run']}")
    files = [item for item in results[0].get("files", []) if item.get("type") == "sra"]
    if not files:
        raise ValueError(f"No full SRA object returned for {row['Run']}")
    file_info = files[0]
    locations = [
        item for item in file_info.get("locations", []) if item.get("service") == "s3"
    ]
    if not locations:
        raise ValueError(f"No public S3 location returned for {row['Run']}")
    url = str(locations[0]["link"]).replace(
        "https://sra-pub-run-odp.s3.amazonaws.com/",
        "https://s3.amazonaws.com/sra-pub-run-odp/",
    )
    return SraFile(
        run_accession=row["Run"],
        biosample=row["BioSample"],
        sample_name=row["SampleName"],
        size_bytes=int(file_info["size"]),
        md5=str(file_info["md5"]).lower(),
        url=url,
    )


def build_sra_manifest(
    rows: list[dict[str, str]], metadata_dir: Path
) -> list[SraFile]:
    """Resolve all SRA objects concurrently and save a reproducible manifest."""
    manifest_path = metadata_dir / "sra_download_manifest.tsv"
    if manifest_path.exists():
        with manifest_path.open(encoding="utf-8", newline="") as handle:
            records = list(csv.DictReader(handle, delimiter="\t"))
        if len(records) != len(rows):
            raise ValueError("Existing SRA manifest row count does not match RunInfo")
        return [
            SraFile(
                run_accession=row["run_accession"],
                biosample=row["biosample"],
                sample_name=row["sample_name"],
                size_bytes=int(row["size_bytes"]),
                md5=row["md5"],
                url=row["url"].replace(
                    "https://sra-pub-run-odp.s3.amazonaws.com/",
                    "https://s3.amazonaws.com/sra-pub-run-odp/",
                ),
            )
            for row in records
        ]
    with ThreadPoolExecutor(max_workers=8) as executor:
        manifest = list(executor.map(locate_sra, rows))
    manifest.sort(key=lambda item: item.run_accession)
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SraFile.__dataclass_fields__)
        for item in manifest:
            writer.writerow(
                (
                    item.run_accession,
                    item.biosample,
                    item.sample_name,
                    item.size_bytes,
                    item.md5,
                    item.url,
                )
            )
    return manifest


def geo_prefix(accession: str) -> str:
    """Return the GEO FTP grouping directory for a Series accession."""
    return f"{accession[:-3]}nnn"


def list_http_files(url: str) -> list[str]:
    """List regular files exposed by an Apache-style directory page."""
    parser = LinkParser()
    parser.feed(run_curl(url, capture=True))
    return sorted(
        {
            link
            for link in parser.links
            if not link.startswith(("/", "?", "#")) and not link.endswith("/")
        }
    )


def prepare_geo_files(accession: str, metadata_dir: Path, geo_dir: Path) -> list[RemoteFile]:
    """Prepare every public GEO family, matrix, and supplementary file."""
    base = f"{GEO_FTP}/{geo_prefix(accession)}/{accession}"
    metadata_urls = {
        f"{accession}_family.soft.gz": f"{base}/soft/{accession}_family.soft.gz",
        f"{accession}_family.xml.tgz": f"{base}/miniml/{accession}_family.xml.tgz",
    }
    files = [
        RemoteFile(url=url, output_path=metadata_dir / name)
        for name, url in metadata_urls.items()
    ]
    for subdirectory in ("matrix", "suppl"):
        directory_url = f"{base}/{subdirectory}/"
        for name in list_http_files(directory_url):
            output_path = (
                metadata_dir / f"{accession}_filelist.txt"
                if name == "filelist.txt"
                else geo_dir / name
            )
            files.append(
                RemoteFile(
                    url=urljoin(directory_url, name),
                    output_path=output_path,
                )
            )
    filelist = next(
        (
            item
            for item in files
            if item.output_path.name == f"{accession}_filelist.txt"
        ),
        None,
    )
    if filelist is not None and not filelist.output_path.exists():
        download_file(filelist)
    if filelist is not None:
        with filelist.output_path.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("#Archive/File") == "Archive":
                    for index, item in enumerate(files):
                        if item.output_path.name == row["Name"]:
                            files[index] = RemoteFile(
                                url=item.url,
                                output_path=item.output_path,
                                size_bytes=int(row["Size"]),
                            )
    return files


def download_many(files: list[RemoteFile], workers: int) -> None:
    """Download independent files with bounded concurrency."""
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(download_file, item): item for item in files}
        for future in as_completed(futures):
            future.result()


def configure_logging(project_dir: Path) -> None:
    """Log key download events to both the terminal and a project log."""
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s")
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    file_handler = logging.FileHandler(
        project_dir / "download.log", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    root_logger.handlers[:] = [file_handler, stream_handler]


def main() -> None:
    """Validate inputs, build manifests, download all files, and verify outputs."""
    args = parse_args()
    accession = args.bioproject.upper()
    if not BIOPROJECT_PATTERN.fullmatch(accession):
        raise ValueError("BioProject accession must match PRJNA followed by digits")
    if not 1 <= args.workers <= 16:
        raise ValueError("--workers must be between 1 and 16")
    if shutil.which("curl.exe" if sys.platform == "win32" else "curl") is None:
        raise RuntimeError("curl is required but was not found on PATH")

    project_dir = args.output_root.resolve() / accession
    metadata_dir = project_dir / "metadata"
    sra_dir = project_dir / "sra"
    geo_dir = project_dir / "geo"
    for directory in (project_dir, metadata_dir, sra_dir, geo_dir):
        directory.mkdir(parents=True, exist_ok=True)
    configure_logging(project_dir)

    logging.info("Resolving BioProject: %s", accession)
    _, geo_accessions = fetch_bioproject_xml(accession, metadata_dir)
    runinfo = fetch_runinfo(accession, metadata_dir)
    manifest = build_sra_manifest(runinfo, metadata_dir)
    total_sra_bytes = sum(item.size_bytes for item in manifest)
    logging.info(
        "SRA target: %d runs, %.2f GiB", len(manifest), total_sra_bytes / 1024**3
    )

    sra_files = [
        RemoteFile(
            url=item.url,
            output_path=sra_dir / f"{item.run_accession}.sra",
            size_bytes=item.size_bytes,
            md5=item.md5,
        )
        for item in manifest
    ]
    download_many(sra_files, args.workers)

    for geo_accession in geo_accessions:
        logging.info("Preparing linked GEO Series: %s", geo_accession)
        geo_files = prepare_geo_files(geo_accession, metadata_dir, geo_dir)
        download_many(geo_files, args.workers)

    logging.info(
        "COMPLETE: %s (%d SRA runs; GEO: %s)",
        accession,
        len(manifest),
        ", ".join(geo_accessions) if geo_accessions else "none linked",
    )


if __name__ == "__main__":
    main()
