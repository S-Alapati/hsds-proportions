"""Thin wrappers around the command line tools the pipeline depends on."""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger(__name__)

REQUIRED = {
    "samtools": "https://www.htslib.org/",
    "filtlong": "https://github.com/rrwick/Filtlong",
    "NanoStat": "https://github.com/wdecoster/nanostat",
}


class ExternalToolError(RuntimeError):
    """Raised when a dependency is missing or exits non-zero."""


def check_tools(names: list[str]) -> list[str]:
    """Return the names that are not on PATH."""
    return [name for name in names if shutil.which(name) is None]


def run(command: list[str], stdout=None) -> subprocess.CompletedProcess:
    log.debug("running: %s", " ".join(str(c) for c in command))
    try:
        return subprocess.run(
            [str(c) for c in command],
            stdout=stdout,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as error:
        message = (error.stderr or "").strip().splitlines()
        tail = message[-1] if message else f"exit status {error.returncode}"
        raise ExternalToolError(f"{command[0]} failed: {tail}") from error


def count_reads(bam: Path, threads: int = 1) -> int:
    result = run(["samtools", "view", "-@", threads, "-c", bam], stdout=subprocess.PIPE)
    return int(result.stdout.strip())


def bam_to_fastq(bam: Path, fastq: Path, threads: int = 1) -> None:
    with open(fastq, "w", encoding="utf-8") as handle:
        run(["samtools", "fastq", "-@", threads, bam], stdout=handle)


def subset_bam(bam: Path, read_ids: Path, output: Path, threads: int = 1) -> None:
    run(["samtools", "view", "-@", threads, "-N", read_ids, "-o", output, bam])


def filtlong_min_quality(fastq: Path, output: Path, min_mean_q: float) -> None:
    with open(output, "w", encoding="utf-8") as handle:
        run(["filtlong", "--min_mean_q", min_mean_q, fastq], stdout=handle)


def read_ids_from_fastq(fastq: Path) -> list[str]:
    """Pull read identifiers from a FASTQ header line.

    Every fourth line is a header, which is safer than grepping for a leading
    ``@`` because that character is also a legal quality score.
    """
    ids: list[str] = []
    with open(fastq, encoding="utf-8") as handle:
        for index, line in enumerate(handle):
            if index % 4 == 0 and line.startswith("@"):
                ids.append(line[1:].split()[0])
    return ids


def nanostat_metric(fastq: Path, metric: str) -> float | None:
    """Return one value from a NanoStat report, or None if it is absent."""
    result = run(["NanoStat", "--fastq", fastq], stdout=subprocess.PIPE)
    for line in result.stdout.splitlines():
        if metric in line:
            return float(line.split()[-1].replace(",", ""))
    return None
