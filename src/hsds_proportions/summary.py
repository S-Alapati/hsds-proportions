"""Allele proportions and the confidence attached to them."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ErrorModel:
    """Combined probability that an allele assignment is correct.

    Three independent terms are multiplied: correct demultiplexing, correct
    basecalling across the informative bases of the motif, and a correct
    methylation call. The result is reported as its complement, an error rate
    per assigned read.

    ``core_bases`` is the number of positions treated as informative. The
    default of five counts the fixed bases either side of the methylated
    adenine rather than the full motif length, on the basis that the spacer is
    unconstrained and a miscall there does not change which family a read
    matches. Set it to whatever suits your motif.
    """

    barcode_accuracy: float = 0.999
    core_bases: int = 5

    def error_rate(self, mean_accuracy_pct: float, min_probability: float) -> float:
        """Per-read error rate as a percentage."""
        base_accuracy = mean_accuracy_pct / 100.0
        methylation_confidence = min_probability / 255.0
        correct = self.barcode_accuracy * (base_accuracy**self.core_bases) * methylation_confidence
        return (1.0 - correct) * 100.0


def proportions(base_counts: dict[str, int]) -> dict[str, float]:
    """Share of assigned bases carried by each allele, as a percentage."""
    total = sum(base_counts.values())
    if total == 0:
        return {name: 0.0 for name in base_counts}
    return {name: count / total * 100.0 for name, count in base_counts.items()}


def format_table(
    base_counts: dict[str, int],
    read_counts: dict[str, int],
    error_rate_pct: float,
    genome_size: int | None = None,
) -> list[str]:
    """Build the summary table as a list of tab separated lines.

    When ``genome_size`` is given, an estimated depth column is added: the
    bases assigned to an allele divided by the genome length. It is a depth
    over the whole genome rather than over the hsdS locus, so it understates
    locus coverage and should be reported as such.
    """
    shares = proportions(base_counts)
    read_total = sum(read_counts.get(name, 0) for name in base_counts)
    header = ["Allele", "Reads", "Total_Bases"]
    if genome_size:
        header.append("Est_Coverage(X)")
    header += ["Read_Proportion(%)", "Base_Proportion(%)", "Error_Rate(%)"]
    lines = ["\t".join(header)]

    for name in sorted(base_counts, key=lambda n: -base_counts[n]):
        reads = read_counts.get(name, 0)
        read_share = reads / read_total * 100.0 if read_total else 0.0
        row = [name, str(reads), str(base_counts[name])]
        if genome_size:
            row.append(f"{base_counts[name] / genome_size:.2f}")
        row += [f"{read_share:.2f}", f"{shares[name]:.2f}", f"{error_rate_pct:.4f}"]
        lines.append("\t".join(row))
    return lines
