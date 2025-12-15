# hsds-proportions

Type I restriction-modification systems with a shufflon-type *hsdS* locus switch
specificity by recombining their target recognition domains, so a clonal culture
is really a mixture of cells carrying different alleles. Each allele leaves a
different 6mA recognition motif on the DNA, which means a long read carrying
methylation calls reports the allele of the cell it came from.

This tool reads a modified-basecalled BAM, assigns each read to an allele from
the motifs it carries, and reports what proportion of the population each allele
accounts for. It was written for *Porphyromonas gingivalis* WW2842, whose *hsdS*
locus produces four alleles (N1C1, N1C2, N2C1 and N2C2), but the motif
definitions are configurable and any comparable locus can be scored.

## What it does

1. Filters reads on mean accuracy with Filtlong, then subsets the original BAM
   so the MM and ML tags survive. FASTQ cannot carry them, which is why the
   filter runs as a round trip rather than in place.
2. Walks the MM and ML tags of every read, keeps the 6mA calls above a
   probability threshold and outside the excluded sequence contexts, and looks
   for the recognition motif of each allele with a methylated adenine at the
   expected offset. Reads are written to one BAM per allele.
3. Totals the bases behind each allele with NanoStat and writes a summary table
   of proportions.

## Installation

The Python package:

```bash
pip install git+https://github.com/S-Alapati/hsds-proportions.git
```

Samtools, Filtlong and NanoStat have to be on PATH. The easiest way to get all
four is conda:

```bash
conda env create -f environment.yml
conda activate hsds-proportions
pip install -e .
```

## Usage

Score every BAM in a directory with the defaults used for WW2842:

```bash
hsds-proportions /path/to/bams -o results
```

Score two files, keep only calls at probability 240 or above, and give samtools
eight threads:

```bash
hsds-proportions barcode05.bam barcode10.bam -o results -p 240 -t 8
```

Process eight barcodes at once, four samtools threads each:

```bash
hsds-proportions /path/to/bams -o results -j 8 -t 4
```

Rescore a BAM that has already been quality filtered:

```bash
hsds-proportions filtered.bam --skip-quality-filter -o results
```

`hsds-proportions --help` lists every option. The ones worth knowing about:

| Option | Default | What it controls |
| --- | --- | --- |
| `-q, --min-quality` | 99.0 | Mean read accuracy as a percentage, passed to Filtlong |
| `-p, --min-probability` | 255 | Lowest 6mA probability accepted, on the 0 to 255 ML scale |
| `--min-spacer`, `--max-spacer` | 7, 7 | Length of the unconstrained run between the two half sites |
| `--motifs` | built in | JSON file of motif families, for a different locus |
| `--exclude` | see below | Sequence context whose calls are discarded, repeatable |
| `--dominance-high/low/switch` | 0.90, 0.75, 5 | How dominant one family has to be before a read is assigned |
| `--genome-size` | 2300000 | Genome length behind the estimated depth column |
| `-t, --threads` | 1 | Threads handed to samtools |
| `-j, --jobs` | 1 | Input files processed at once |

## Output

One directory per input BAM, named after the sample and the thresholds used, for
example `barcode05_Q99_P255/`:

| File | Contents |
| --- | --- |
| `*_Q99_filtered.bam` | Reads that passed the accuracy filter |
| `*_N1C1.bam` and the other three | Reads assigned to each allele |
| `*_Ambiguous.bam` | Reads carrying motifs from more than one family with no clear majority |
| `*_Other.bam` | Reads with no 6mA calls, no motifs, or no MM tag |
| `*_motifs.tsv` | Every motif occurrence, with its sequence and methylation pattern |
| `*_summary_report.tsv` | Read and base counts per allele, proportions, error rate |

In the methylation pattern, `M` marks a retained call, `m` a call dropped for
sitting in an excluded context, and `.` an unmodified base.

## How a read is assigned

A read carrying motifs from one family only goes to that family. Where more than
one family is represented, one of them has to account for most of the
occurrences: 90% on reads with at least five motifs, 75% below that, on the
grounds that a handful of motifs is weak evidence and a single miscalled base
should not decide the allele. Reads that clear neither bar are written to the
ambiguous file rather than forced into a call, and they are left out of the
proportions.

## Scoring a different locus

Motif families are a JSON mapping, one entry per allele, with the spacer written
as `{spacer}` and the offsets at which the methylated adenine is expected:

```json
{
  "N1C1": {"fwd": "TCA{spacer}TGT", "rev": "ACA{spacer}TGA", "offsets": [2]},
  "N2C2": {"fwd": "GTA{spacer}TTA", "rev": "TAA{spacer}TAC", "offsets": [2]}
}
```

Pass it with `--motifs my_motifs.json`. `examples/ww2842_motifs.json` holds the
built in definitions as a starting point.

Three sequence contexts are excluded by default, `CGCAG` with the adenine at
offset 3, `CCAGG` at offset 2 and `GGACC` at offset 2. They produce 6mA calls
that have nothing to do with the Type I system and would otherwise be counted.
Add others with `--exclude SEQ:OFFSET`, or clear the list with
`--no-default-exclusions`.

## The error rate column

The summary reports a single error rate per assigned read, the complement of
three probabilities multiplied together: correct demultiplexing, correct
basecalling across the informative bases of the motif, and a correct methylation
call. With the defaults, `0.999 x 0.99^5 x 1.0`, this gives 4.9961%.

Two things to be aware of before quoting it. The exponent is the number of
fixed bases flanking the methylated adenine rather than the full motif length,
on the basis that a miscall in the spacer does not change which family a read
matches. And the methylation term is the probability threshold expressed as a
fraction of 255, which at the default of 255 contributes nothing. Both are set
by `--error-model-bases` and `--barcode-accuracy` if a different treatment suits
your data better.

The estimated depth column divides the bases assigned to an allele by the whole
genome length, so it is a genome-equivalent depth and understates coverage at
the locus itself. Set `--genome-size 0` to leave it out.

## Known behaviour

Motif matches within a pattern are non-overlapping, which is the behaviour of
`re.finditer`. Two overlapping copies of the same motif count once.

By default the methylation pattern marks every raw 6mA call, including those
below the probability threshold, which is how the original analysis was run.
`--annotate-kept-only` restricts the marks to calls that passed. This affects
the pattern string in the TSV only, never which reads are assigned where.

Reads are assigned on the motifs they happen to span, so short reads carry less
evidence than long ones. Proportions are reported on a base basis as well as a
read basis for that reason, and the two are worth comparing.

## Requirements

Python 3.9 or later with pysam, plus samtools, Filtlong and NanoStat on PATH.

## Citation

If this is useful in published work please cite the repository and the
manuscript it was written for. See `CITATION.cff`.

## Licence

MIT. See `LICENSE`.
