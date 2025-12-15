"""Running the three stages over one BAM file."""

from __future__ import annotations

import logging
from contextlib import ExitStack
from dataclasses import dataclass, field
from pathlib import Path

from . import external
from .assign import AMBIGUOUS, UNASSIGNED, assign_allele, find_motifs
from .modifications import annotation_string, filter_calls, parse_mm_ml
from .motifs import MotifSet
from .summary import ErrorModel, format_table

log = logging.getLogger(__name__)


@dataclass
class Settings:
    """Everything the pipeline needs that is not the input file itself."""

    motif_set: MotifSet
    output_dir: Path
    mod_code: str = "A+a."
    min_probability: float = 255.0
    min_mean_quality: float = 99.0
    threads: int = 1
    genome_size: int | None = 2_300_000
    error_model: ErrorModel = field(default_factory=ErrorModel)
    high_threshold: float = 0.90
    low_threshold: float = 0.75
    switch_at: int = 5
    skip_quality_filter: bool = False
    write_summary: bool = True
    keep_intermediate: bool = False
    annotate_kept_only: bool = False
    flat_output: bool = False


@dataclass
class SampleResult:
    sample: str
    output_dir: Path
    reads_total: int = 0
    reads_passed: int = 0
    read_counts: dict[str, int] = field(default_factory=dict)
    base_counts: dict[str, int] = field(default_factory=dict)
    summary_path: Path | None = None
    failed: str | None = None


def quality_filter(bam: Path, settings: Settings, prefix: Path) -> tuple[Path, int, int]:
    """Keep reads whose mean accuracy clears the threshold.

    Filtlong works on FASTQ, so the BAM is converted, filtered, and the
    surviving identifiers are used to subset the original BAM. Going back to
    the BAM matters because the FASTQ has lost the MM and ML tags.
    """
    filtered_bam = prefix.with_name(f"{prefix.name}_Q{settings.min_mean_quality:g}_filtered.bam")
    fastq = prefix.with_name(f"{prefix.name}.tmp.fastq")
    fastq_passed = prefix.with_name(f"{prefix.name}.tmp.pass.fastq")
    id_list = prefix.with_name(f"{prefix.name}.tmp.ids.txt")

    try:
        total = external.count_reads(bam, settings.threads)
        log.info("converting %s to FASTQ", bam.name)
        external.bam_to_fastq(bam, fastq, settings.threads)

        log.info("filtering at mean accuracy >= %g%%", settings.min_mean_quality)
        external.filtlong_min_quality(fastq, fastq_passed, settings.min_mean_quality)

        ids = external.read_ids_from_fastq(fastq_passed)
        id_list.write_text("\n".join(ids), encoding="utf-8")
        log.info("%d of %d reads passed, %d discarded", len(ids), total, total - len(ids))

        if not ids:
            raise external.ExternalToolError("no reads passed the quality filter")

        external.subset_bam(bam, id_list, filtered_bam, settings.threads)
        return filtered_bam, total, len(ids)
    finally:
        if not settings.keep_intermediate:
            for path in (fastq, fastq_passed, id_list):
                path.unlink(missing_ok=True)


def split_by_allele(bam: Path, settings: Settings, prefix: Path) -> tuple[dict[str, Path], dict[str, int]]:
    """Score every read and write it to the file for its allele."""
    import pysam  # imported here so the pure logic stays importable without it

    families = list(settings.motif_set.family_names)
    categories = families + [AMBIGUOUS, UNASSIGNED]
    paths = {name: prefix.with_name(f"{prefix.name}_{name}.bam") for name in categories}
    counts = {name: 0 for name in categories}
    motif_tsv = prefix.with_name(f"{prefix.name}_motifs.tsv")

    with ExitStack() as stack:
        source = stack.enter_context(
            pysam.AlignmentFile(str(bam), "rb", check_sq=False, require_index=False)
        )
        handles = {
            name: stack.enter_context(pysam.AlignmentFile(str(path), "wb", template=source))
            for name, path in paths.items()
        }
        tsv = stack.enter_context(open(motif_tsv, "w", encoding="utf-8"))
        tsv.write("Read_ID\tMotif\tSequence\tMethylation\tAssigned_To\n")

        for read in source:
            sequence = read.query_sequence
            if not sequence:
                handles[UNASSIGNED].write(read)
                counts[UNASSIGNED] += 1
                continue

            try:
                mm_tag = read.get_tag("MM")
                ml_tag = read.get_tag("ML")
            except KeyError:
                handles[UNASSIGNED].write(read)
                counts[UNASSIGNED] += 1
                continue

            if settings.mod_code not in mm_tag:
                handles[UNASSIGNED].write(read)
                counts[UNASSIGNED] += 1
                continue

            calls = parse_mm_ml(sequence, mm_tag, ml_tag)
            positions = filter_calls(
                sequence, calls, settings.mod_code, settings.min_probability, settings.motif_set.exclusions
            )
            if not positions.kept:
                handles[UNASSIGNED].write(read)
                counts[UNASSIGNED] += 1
                continue

            annotation = annotation_string(len(sequence), positions, settings.annotate_kept_only)
            hits = find_motifs(sequence, positions.kept, annotation, settings.motif_set)
            allele = assign_allele(
                hits,
                settings.motif_set.family_names,
                settings.high_threshold,
                settings.low_threshold,
                settings.switch_at,
            )

            handles[allele].write(read)
            counts[allele] += 1
            if allele != UNASSIGNED:
                for hit in hits:
                    tsv.write(
                        f"{read.query_name}\t{hit.pattern}\t{hit.sequence}\t"
                        f"{hit.methylation}\t{allele}\n"
                    )

    for name, count in counts.items():
        log.info("%-10s %7d reads -> %s", name, count, paths[name].name)
    return {name: paths[name] for name in families}, counts


def base_totals(allele_bams: dict[str, Path], settings: Settings, prefix: Path) -> dict[str, int]:
    """Total sequenced bases behind each allele, via NanoStat."""
    totals: dict[str, int] = {}
    for name, bam in allele_bams.items():
        if not bam.exists() or bam.stat().st_size == 0:
            totals[name] = 0
            continue
        fastq = prefix.with_name(f"{prefix.name}.tmp.{name}.fastq")
        try:
            external.bam_to_fastq(bam, fastq, settings.threads)
            value = external.nanostat_metric(fastq, "Total bases:")
            totals[name] = int(value) if value else 0
        except external.ExternalToolError as error:
            log.warning("NanoStat failed on %s: %s", bam.name, error)
            totals[name] = 0
        finally:
            fastq.unlink(missing_ok=True)
    return totals


def process_bam(bam: Path, settings: Settings) -> SampleResult:
    """Run the whole pipeline for one input BAM."""
    sample = bam.stem
    if settings.flat_output:
        out_dir = settings.output_dir
    else:
        tag = f"Q{settings.min_mean_quality:g}_P{settings.min_probability:g}"
        out_dir = settings.output_dir / f"{sample}_{tag}"
    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / sample
    result = SampleResult(sample=sample, output_dir=out_dir)
    log.info("writing results for %s to %s", sample, out_dir)

    try:
        if settings.skip_quality_filter:
            log.info("skipping the quality filter, scoring %s as given", bam.name)
            working_bam = bam
            result.reads_total = result.reads_passed = external.count_reads(bam, settings.threads)
        else:
            working_bam, result.reads_total, result.reads_passed = quality_filter(bam, settings, prefix)

        allele_bams, result.read_counts = split_by_allele(working_bam, settings, prefix)

        if settings.write_summary:
            result.base_counts = base_totals(allele_bams, settings, prefix)
            result.summary_path = write_summary(result, settings, prefix)
    except (external.ExternalToolError, OSError, ValueError) as error:
        log.error("%s: %s", sample, error)
        result.failed = str(error)

    return result


def write_summary(result: SampleResult, settings: Settings, prefix: Path) -> Path:
    """Write the per-sample proportions table."""
    path = prefix.with_name(f"{prefix.name}_summary_report.tsv")
    error_rate = settings.error_model.error_rate(settings.min_mean_quality, settings.min_probability)
    assigned_reads = sum(result.read_counts.get(n, 0) for n in result.base_counts)
    table = format_table(result.base_counts, result.read_counts, error_rate, settings.genome_size)

    with open(path, "w", encoding="utf-8") as handle:
        handle.write(f"Total_Reads_Analysed: {result.reads_total}\n")
        handle.write(f"Reads_Passing_Quality_Filter: {result.reads_passed}\n")
        handle.write(f"Reads_Assigned_To_An_Allele: {assigned_reads}\n")
        handle.write(f"Total_Bases_Assigned: {sum(result.base_counts.values())}\n")
        handle.write(
            f"Parameters: min_mean_quality={settings.min_mean_quality:g} "
            f"min_probability={settings.min_probability:g} "
            f"spacer={settings.motif_set.min_spacer}-{settings.motif_set.max_spacer}\n\n"
        )
        handle.write("\n".join(table) + "\n")

    print("\n".join(table))
    log.info("summary written to %s", path)
    return path
