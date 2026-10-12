#!/usr/bin/env python3
"""Query NCBI and convert BioProject SRA files to sample-title FASTQ on Linux."""

from __future__ import annotations

import argparse
import csv
import io
import json
import logging
import re
import shutil
import subprocess
import tempfile
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import quote
from urllib.request import Request, urlopen


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BIOPROJECT_PATTERN = re.compile(r"^PRJNA[0-9]+$")
SAFE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9._-]+$")


@dataclass(frozen=True)
class Sample:
    """Store a run, its NCBI identifiers, and biological sample title."""

    run_accession: str
    sample_name: str
    sample_accession: str
    sample_title: str

    @property
    def output_prefix(self) -> str:
        """Return the biological sample title used as the FASTQ prefix."""
        return self.sample_title


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Query NCBI and convert paired SRA to sample-title FASTQ.gz."
    )
    parser.add_argument("bioproject", help="NCBI accession, e.g. PRJNA1147209")
    parser.add_argument("--sra-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--temp-dir", type=Path, required=True)
    parser.add_argument("--metadata-dir", type=Path)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--threads-per-job", type=int, default=16)
    return parser.parse_args()


def run_command(command: list[str]) -> None:
    """Run one external command and raise on failure."""
    logging.info("RUN: %s", " ".join(command))
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {' '.join(command)}"
        )


def check_programs() -> None:
    """Require conversion programs before querying or writing output."""
    missing = [
        name
        for name in ("vdb-validate", "fasterq-dump", "pigz")
        if shutil.which(name) is None
    ]
    if missing:
        raise RuntimeError(f"Missing required programs: {', '.join(missing)}")


def fetch_text(url: str, attempts: int = 3) -> str:
    """Retrieve verified HTTPS text, retrying transport failures only."""
    request = Request(url, headers={"User-Agent": "ConvertSraToFastq/1.0"})
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(request, timeout=60) as response:
                return response.read().decode("utf-8")
        except Exception as error:
            last_error = error
            if attempt < attempts:
                time.sleep(5 * attempt)
    raise RuntimeError(
        f"Failed to retrieve NCBI URL after {attempts} attempts: {url}"
    ) from last_error


def batched(values: list[str], size: int) -> Iterable[list[str]]:
    """Yield fixed-size batches."""
    for start in range(0, len(values), size):
        yield values[start : start + size]


def query_sample_titles(sra_ids: list[str]) -> dict[str, str]:
    """Map run accessions to biological titles from SRA Experiment XML."""
    titles: dict[str, str] = {}
    for id_batch in batched(sra_ids, 200):
        url = f"{EUTILS}/esummary.fcgi?db=sra&id={','.join(id_batch)}&retmode=json"
        data: dict[str, Any] = json.loads(fetch_text(url))
        result = data["result"]
        for uid in result["uids"]:
            record = result[uid]
            exp_root = ET.fromstring(f"<Root>{record['expxml']}</Root>")
            experiment = exp_root.find("Experiment")
            experiment_name = "" if experiment is None else experiment.get("name", "")
            if ":" not in experiment_name or ";" not in experiment_name:
                raise ValueError(f"Cannot parse sample title from: {experiment_name}")
            title = experiment_name.split(":", 1)[1].split(";", 1)[0].strip()
            if not title:
                raise ValueError(f"Empty sample title in: {experiment_name}")
            runs_root = ET.fromstring(f"<Root>{record['runs']}</Root>")
            run_accessions = [node.get("acc", "").strip() for node in runs_root.findall("Run")]
            if not run_accessions or any(not value for value in run_accessions):
                raise ValueError(f"No valid run accession in SRA summary {uid}")
            for run_accession in run_accessions:
                if run_accession in titles:
                    raise ValueError(f"Duplicate title mapping: {run_accession}")
                titles[run_accession] = title
    return titles


def query_sra_runinfo(bioproject: str) -> tuple[list[Sample], str]:
    """Query RunInfo plus SRA titles and return normalized metadata."""
    search_url = (
        f"{EUTILS}/esearch.fcgi?db=sra"
        f"&term={quote(bioproject + '[BioProject]')}&retmode=json&retmax=100000"
    )
    search: dict[str, Any] = json.loads(fetch_text(search_url))
    result = search["esearchresult"]
    ids = [str(value) for value in result["idlist"]]
    expected_count = int(result["count"])
    if not ids or len(ids) != expected_count:
        raise ValueError(
            f"Expected {expected_count} SRA records, retrieved {len(ids)}"
        )

    rows: list[dict[str, str]] = []
    required = {"Run", "SampleName", "Sample", "LibraryLayout", "BioProject"}
    for id_batch in batched(ids, 200):
        url = (
            f"{EUTILS}/efetch.fcgi?db=sra&id={','.join(id_batch)}"
            "&rettype=runinfo&retmode=text"
        )
        reader = csv.DictReader(io.StringIO(fetch_text(url)))
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("NCBI RunInfo is missing required columns")
        rows.extend(dict(row) for row in reader)
    if len(rows) != expected_count:
        raise ValueError(f"Expected {expected_count} RunInfo rows, received {len(rows)}")
    if {row["BioProject"] for row in rows} != {bioproject}:
        raise ValueError("RunInfo contains unexpected BioProject accessions")
    non_paired = [row["Run"] for row in rows if row["LibraryLayout"] != "PAIRED"]
    if non_paired:
        raise ValueError(f"Non-paired runs found: {', '.join(non_paired)}")

    titles = query_sample_titles(ids)
    run_accessions = {row["Run"].strip() for row in rows}
    if set(titles) != run_accessions:
        raise ValueError("SRA title mappings do not match RunInfo accessions")
    samples = [
        Sample(
            run_accession=row["Run"].strip(),
            sample_name=row["SampleName"].strip(),
            sample_accession=row["Sample"].strip(),
            sample_title=titles[row["Run"].strip()],
        )
        for row in rows
    ]
    validate_sample_mappings(samples)
    samples.sort(key=lambda item: item.run_accession)

    output = io.StringIO()
    writer = csv.writer(output, delimiter="\t", lineterminator="\n")
    writer.writerow(("run_accession", "sample_name", "sample", "sample_title"))
    for sample in samples:
        writer.writerow(
            (
                sample.run_accession,
                sample.sample_name,
                sample.sample_accession,
                sample.sample_title,
            )
        )
    return samples, output.getvalue()


def validate_sample_mappings(samples: list[Sample]) -> None:
    """Validate identifiers, filename safety, and title uniqueness."""
    if not samples:
        raise ValueError("NCBI returned no samples")
    for sample in samples:
        values = (
            sample.run_accession,
            sample.sample_name,
            sample.sample_accession,
            sample.sample_title,
        )
        if any(not value for value in values):
            raise ValueError(f"Empty identifier in mapping: {sample}")
        if any(SAFE_NAME_PATTERN.fullmatch(value) is None for value in values):
            raise ValueError(f"Unsafe identifier for filename: {sample}")
    if len({item.run_accession for item in samples}) != len(samples):
        raise ValueError("Duplicate run accessions")
    if len({item.sample_title for item in samples}) != len(samples):
        raise ValueError("Biological sample titles are not unique")


def save_metadata(path: Path, content: str) -> None:
    """Save metadata without replacing different existing content."""
    if path.exists():
        if path.read_text(encoding="utf-8") != content:
            raise ValueError(f"Existing metadata differs from NCBI query: {path}")
        logging.info("NCBI metadata is unchanged: %s", path)
        return
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write(content)


def validate_paths(samples: list[Sample], sra_dir: Path, output_dir: Path) -> None:
    """Check all input and final paths before conversion."""
    missing = [
        str(sra_dir / f"{item.run_accession}.sra")
        for item in samples
        if not (sra_dir / f"{item.run_accession}.sra").is_file()
        or (sra_dir / f"{item.run_accession}.sra").stat().st_size == 0
    ]
    if missing:
        raise ValueError(f"Missing or empty SRA files ({len(missing)}):\n" + "\n".join(missing[:10]))
    existing = [
        str(path)
        for item in samples
        for path in (
            output_dir / f"{item.output_prefix}_R1.fastq.gz",
            output_dir / f"{item.output_prefix}_R2.fastq.gz",
        )
        if path.exists()
    ]
    if existing:
        raise FileExistsError("Refusing to overwrite FASTQ files:\n" + "\n".join(existing[:10]))


def convert_sample(
    sample: Sample, sra_dir: Path, output_dir: Path, temp_root: Path, threads: int
) -> None:
    """Validate, convert, compress, test, and rename one paired SRA run."""
    sra_path = sra_dir / f"{sample.run_accession}.sra"
    final_r1 = output_dir / f"{sample.output_prefix}_R1.fastq.gz"
    final_r2 = output_dir / f"{sample.output_prefix}_R2.fastq.gz"
    run_command(["vdb-validate", str(sra_path)])
    with tempfile.TemporaryDirectory(prefix=f"{sample.run_accession}_", dir=temp_root) as tmp:
        work_dir = Path(tmp)
        tool_tmp = work_dir / "fasterq_tmp"
        tool_tmp.mkdir()
        run_command(
            [
                "fasterq-dump", "--split-files", "--threads", str(threads),
                "--temp", str(tool_tmp), "--outdir", str(work_dir), str(sra_path),
            ]
        )
        fastq_r1 = work_dir / f"{sample.run_accession}_1.fastq"
        fastq_r2 = work_dir / f"{sample.run_accession}_2.fastq"
        if any(not path.is_file() or path.stat().st_size == 0 for path in (fastq_r1, fastq_r2)):
            raise ValueError(f"Missing or empty paired FASTQ for {sample.run_accession}")
        run_command(["pigz", "-p", str(threads), str(fastq_r1), str(fastq_r2)])
        gzip_r1 = Path(f"{fastq_r1}.gz")
        gzip_r2 = Path(f"{fastq_r2}.gz")
        run_command(["pigz", "-t", str(gzip_r1)])
        run_command(["pigz", "-t", str(gzip_r2)])
        gzip_r1.rename(final_r1)
        gzip_r2.rename(final_r2)


def main() -> None:
    """Query NCBI, validate inputs, and convert all runs concurrently."""
    args = parse_args()
    bioproject = args.bioproject.upper()
    if BIOPROJECT_PATTERN.fullmatch(bioproject) is None:
        raise ValueError("BioProject must match PRJNA followed by digits")
    if args.jobs < 1 or args.threads_per_job < 1:
        raise ValueError("--jobs and --threads-per-job must be positive")
    check_programs()
    sra_dir = args.sra_dir.resolve()
    output_dir = args.output_dir.resolve()
    temp_dir = args.temp_dir.resolve()
    metadata_dir = args.metadata_dir.resolve() if args.metadata_dir else output_dir / "metadata"
    for directory in (output_dir, temp_dir, metadata_dir):
        directory.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s\t%(levelname)s\t%(message)s",
        handlers=[logging.FileHandler(output_dir / "convert_sra_to_fastq.log"), logging.StreamHandler()],
    )
    samples, metadata = query_sra_runinfo(bioproject)
    save_metadata(metadata_dir / f"{bioproject}_sra_runinfo.tsv", metadata)
    validate_paths(samples, sra_dir, output_dir)
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [
            executor.submit(
                convert_sample, item, sra_dir, output_dir, temp_dir, args.threads_per_job
            )
            for item in samples
        ]
        for future in as_completed(futures):
            future.result()
    logging.info("All %d samples completed", len(samples))


if __name__ == "__main__":
    main()
