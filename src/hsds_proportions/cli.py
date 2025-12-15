"""Command line entry point."""

from __future__ import annotations

import argparse
import logging
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from . import __version__, external
from .motifs import build_motif_set, list_presets, load_definition, load_preset, parse_exclusion
from .pipeline import Settings, process_bam
from .summary import ErrorModel

log = logging.getLogger("hsds")

OUTPUT_SUFFIXES = ("_filtered.bam", "_motifs.tsv", "_summary_report.tsv")


def collect_inputs(paths: list[Path], pattern: str) -> list[Path]:
    """Expand the given files and directories into a list of BAMs.

    Anything this tool has written before is skipped, so rerunning in place
    does not feed the outputs back in as inputs.
    """
    found: list[Path] = []
    for path in paths:
        if path.is_dir():
            found.extend(sorted(path.glob(pattern)))
        elif path.is_file():
            found.append(path)
        else:
            log.warning("%s does not exist, skipping", path)

    keep = []
    for bam in found:
        if bam.name.endswith(OUTPUT_SUFFIXES) or bam.stem.endswith(("_Ambiguous", "_Other")):
            log.debug("skipping %s, it looks like an output of a previous run", bam.name)
            continue
        keep.append(bam)
    return keep


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hsds-proportions",
        description=(
            "Assign long reads to hsdS alleles from the 6mA motifs they carry, "
            "and report the proportion of the population each allele accounts for."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("inputs", nargs="*", type=Path, default=[Path(".")],
                        help="BAM files, or directories to search")
    parser.add_argument("--pattern", default="*.bam",
                        help="glob used when an input is a directory")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("hsds_results"),
                        help="where results are written")
    parser.add_argument("--flat", action="store_true",
                        help="write into the output directory directly instead of one subdirectory per sample")

    thresholds = parser.add_argument_group("thresholds")
    thresholds.add_argument("-q", "--min-quality", type=float, default=99.0, metavar="PCT",
                            help="minimum mean read accuracy passed to filtlong, as a percentage")
    thresholds.add_argument("-p", "--min-probability", type=float, default=255.0, metavar="P",
                            help="minimum 6mA call probability on the 0 to 255 ML scale")
    thresholds.add_argument("--skip-quality-filter", action="store_true",
                            help="score the input as given, for BAMs that are already filtered")

    motif = parser.add_argument_group("motif definition")
    motif.add_argument("--motifs", type=Path, metavar="JSON",
                       help="motif definition for your organism, see --list-presets for a skeleton")
    motif.add_argument("--preset", metavar="NAME",
                       help=f"bundled motif definition, one of: {', '.join(list_presets())}")
    motif.add_argument("--list-presets", action="store_true",
                       help="describe the bundled presets and exit")
    motif.add_argument("--min-spacer", type=int, metavar="N",
                       help="override the shortest spacer for every family")
    motif.add_argument("--max-spacer", type=int, metavar="N",
                       help="override the longest spacer for every family")
    motif.add_argument("--exclude", action="append", metavar="SEQ:OFFSET",
                       help="sequence context whose calls are discarded, repeatable")
    motif.add_argument("--no-default-exclusions", action="store_true",
                       help="drop the exclusion contexts that come with the preset")
    motif.add_argument("--mod-code", default="A+a.",
                       help="modification code to read from the MM tag")
    motif.add_argument("--annotate-kept-only", action="store_true",
                       help="mark only calls above the probability threshold in the methylation string")

    calling = parser.add_argument_group("allele calling")
    calling.add_argument("--dominance-high", type=float, default=0.90, metavar="F",
                         help="share one family needs on reads with many motifs")
    calling.add_argument("--dominance-low", type=float, default=0.75, metavar="F",
                         help="share one family needs on reads with few motifs")
    calling.add_argument("--dominance-switch", type=int, default=5, metavar="N",
                         help="motif count at which the higher share is required")

    reporting = parser.add_argument_group("reporting")
    reporting.add_argument("--genome-size", type=int, metavar="BP",
                           help="genome length for the estimated depth column, 0 to omit it "
                                "(defaults to the value in the motif definition, if it has one)")
    reporting.add_argument("--barcode-accuracy", type=float, default=0.999, metavar="P",
                           help="probability a read is demultiplexed correctly")
    reporting.add_argument("--error-model-bases", type=int, default=5, metavar="N",
                           help="informative bases of the motif used in the error model")
    reporting.add_argument("--no-summary", action="store_true",
                           help="skip the NanoStat step and write no summary table")

    runtime = parser.add_argument_group("runtime")
    runtime.add_argument("-t", "--threads", type=int, default=1,
                         help="threads handed to samtools")
    runtime.add_argument("-j", "--jobs", type=int, default=1,
                         help="input files processed at once")
    runtime.add_argument("--keep-intermediate", action="store_true",
                         help="keep the temporary FASTQ and read identifier files")
    runtime.add_argument("-v", "--verbose", action="store_true", help="show debug messages")
    runtime.add_argument("--quiet", action="store_true", help="warnings and errors only")
    runtime.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return parser


def settings_from_args(args: argparse.Namespace) -> Settings:
    if args.motifs:
        definition = load_definition(args.motifs)
        label = str(args.motifs)
    else:
        definition = load_preset(args.preset)
        label = args.preset

    families = definition["families"]
    if args.min_spacer is not None or args.max_spacer is not None:
        for spec in families.values():
            spec.pop("spacer", None)
            if args.min_spacer is not None:
                spec["min_spacer"] = args.min_spacer
            if args.max_spacer is not None:
                spec["max_spacer"] = args.max_spacer

    exclusions = [] if args.no_default_exclusions else [tuple(e) for e in definition.get("exclusions", [])]
    for text in args.exclude or []:
        exclusions.append(parse_exclusion(text))

    motif_set = build_motif_set(families, exclusions=exclusions, name=label)
    log.info("scoring with %s", motif_set.describe())

    genome_size = definition.get("genome_size") if args.genome_size is None else args.genome_size
    return Settings(
        motif_set=motif_set,
        output_dir=args.output_dir,
        mod_code=args.mod_code,
        min_probability=args.min_probability,
        min_mean_quality=args.min_quality,
        threads=max(1, args.threads),
        genome_size=genome_size or None,
        error_model=ErrorModel(args.barcode_accuracy, args.error_model_bases),
        high_threshold=args.dominance_high,
        low_threshold=args.dominance_low,
        switch_at=args.dominance_switch,
        skip_quality_filter=args.skip_quality_filter,
        write_summary=not args.no_summary,
        keep_intermediate=args.keep_intermediate,
        annotate_kept_only=args.annotate_kept_only,
        flat_output=args.flat,
    )


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    level = logging.DEBUG if args.verbose else logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(level=level, format="%(levelname)s %(message)s")

    if args.list_presets:
        for name in list_presets():
            definition = load_preset(name)
            print(f"{name}\n    {definition.get('name', '')}")
            print(f"    alleles: {', '.join(definition['families'])}")
            if definition.get("reference"):
                print(f"    {definition['reference']}")
        return 0

    if not args.preset and not args.motifs:
        log.error("give a motif definition for your organism with --motifs")
        log.error("run --list-presets for a skeleton to copy, or see examples/")
        return 2

    needed = ["samtools"]
    if not args.skip_quality_filter:
        needed.append("filtlong")
    if not args.no_summary:
        needed.append("NanoStat")
    missing = external.check_tools(needed)
    if missing:
        for name in missing:
            log.error("%s is not on PATH, see %s", name, external.REQUIRED.get(name, ""))
        return 2

    bams = collect_inputs(args.inputs, args.pattern)
    if not bams:
        log.error("no input BAM files found")
        return 1
    log.info("found %d BAM file(s)", len(bams))

    try:
        settings = settings_from_args(args)
    except (ValueError, OSError, KeyError) as error:
        log.error("%s", error)
        return 2
    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.jobs > 1 and len(bams) > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            results = list(pool.map(process_bam, bams, [settings] * len(bams)))
    else:
        results = [process_bam(bam, settings) for bam in bams]

    failures = [r for r in results if r.failed]
    for result in failures:
        log.error("%s did not complete: %s", result.sample, result.failed)
    log.info("%d of %d samples completed", len(results) - len(failures), len(results))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
