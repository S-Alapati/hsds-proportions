"""Finding motifs on a read and deciding which allele the read came from."""

from __future__ import annotations

from dataclasses import dataclass

from .motifs import MotifSet

#: Reads that carry no scorable motif, and reads whose motifs disagree.
UNASSIGNED = "Other"
AMBIGUOUS = "Ambiguous"


@dataclass(frozen=True)
class MotifHit:
    family: str
    pattern: str
    sequence: str
    methylation: str


def find_motifs(
    sequence: str,
    methylated_positions: set[int],
    annotation: str,
    motif_set: MotifSet,
) -> list[MotifHit]:
    """Return every motif occurrence that carries a methylated adenine.

    Matches are non-overlapping within a pattern, which is the behaviour of
    ``re.finditer`` and is deliberate: two overlapping copies of the same motif
    would otherwise be counted twice when deciding the dominant family.
    """
    hits: list[MotifHit] = []
    for label, (pattern, offsets) in motif_set.patterns.items():
        family = label.split("_")[0]
        for match in pattern.finditer(sequence):
            start, end = match.start(), match.end()
            if any(start + offset in methylated_positions for offset in offsets):
                hits.append(
                    MotifHit(
                        family=family,
                        pattern=label,
                        sequence=sequence[start:end],
                        methylation=annotation[start:end],
                    )
                )
    return hits


def assign_allele(
    hits: list[MotifHit],
    families: tuple[str, ...],
    high_threshold: float = 0.90,
    low_threshold: float = 0.75,
    switch_at: int = 5,
) -> str:
    """Decide which allele a read belongs to from the motifs it carries.

    A read with motifs from one family only is assigned to it. Otherwise one
    family has to account for most of the occurrences, and the bar is raised on
    reads with enough motifs to support it: ``high_threshold`` once there are
    at least ``switch_at`` occurrences, ``low_threshold`` below that. Reads that
    clear neither are held as ambiguous rather than forced into a call.
    """
    if not hits:
        return UNASSIGNED

    counts = {family: 0 for family in families}
    for hit in hits:
        if hit.family in counts:
            counts[hit.family] += 1

    present = [family for family, count in counts.items() if count]
    if not present:
        return UNASSIGNED
    if len(present) == 1:
        return present[0]

    total = sum(counts.values())
    threshold = high_threshold if total >= switch_at else low_threshold
    winners = [f for f in present if counts[f] / total >= threshold]
    return winners[0] if len(winners) == 1 else AMBIGUOUS
