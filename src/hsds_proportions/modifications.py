"""Reading 6mA calls out of the SAM MM and ML tags.

The MM tag stores, for each modified base type, the number of unmodified bases
of that type to skip before the next modified one. ML holds the matching
probabilities as bytes on a 0 to 255 scale. Walking the two together gives the
position of every call on the read together with its confidence.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

log = logging.getLogger(__name__)


@dataclass
class ModifiedPositions:
    """6mA calls on a single read, split by how they were treated."""

    kept: set[int]
    excluded: set[int]
    all_calls: set[int]


def parse_mm_ml(read_seq: str, mm_tag: str, ml_tag) -> dict[int, tuple[str, int]]:
    """Map read position to (modification code, probability).

    Where the same position is reported more than once, which happens when a
    basecaller emits several modification codes for one base, the call with the
    highest probability wins.
    """
    try:
        probabilities = list(ml_tag)
    except TypeError:
        log.debug("ML tag is not iterable, skipping read")
        return {}

    calls: dict[int, tuple[str, int]] = {}
    base_positions: dict[str, list[int]] = {}
    consumed = 0

    # A leading "Mm:Z:" is present on some inputs, so take the payload only.
    payload = mm_tag.split(":", 2)[-1]

    for group in payload.split(";"):
        if not group:
            continue
        header, *deltas_text = group.split(",")
        base_type = header[:1]
        if base_type not in "ACGTN":
            continue

        deltas = [int(d) for d in deltas_text if d]
        if consumed + len(deltas) > len(probabilities):
            log.debug("ML tag shorter than MM implies, truncating group %s", header)
            continue
        group_probs = probabilities[consumed : consumed + len(deltas)]
        consumed += len(deltas)

        if base_type not in base_positions:
            base_positions[base_type] = [i for i, b in enumerate(read_seq) if b == base_type]
        positions = base_positions[base_type]

        cursor = 0
        for delta, probability in zip(deltas, group_probs):
            cursor += delta
            if cursor >= len(positions):
                break
            position = positions[cursor]
            if position not in calls or probability > calls[position][1]:
                calls[position] = (header, probability)
            cursor += 1

    return calls


def filter_calls(
    read_seq: str,
    calls: dict[int, tuple[str, int]],
    mod_code: str,
    min_probability: float,
    exclusions: tuple[tuple[str, int], ...],
) -> ModifiedPositions:
    """Drop calls that fall in an excluded context or below the threshold."""
    kept: set[int] = set()
    excluded: set[int] = set()
    seen: set[int] = set()
    length = len(read_seq)

    for position, (code, probability) in calls.items():
        if code != mod_code:
            continue
        seen.add(position)

        in_excluded_context = False
        for context, offset in exclusions:
            start = position - offset
            end = start + len(context)
            if start >= 0 and end <= length and read_seq[start:end] == context:
                in_excluded_context = True
                break

        if in_excluded_context:
            excluded.add(position)
        elif probability >= min_probability:
            kept.add(position)

    return ModifiedPositions(kept=kept, excluded=excluded, all_calls=seen)


def annotation_string(length: int, positions: ModifiedPositions, kept_only: bool = False) -> str:
    """Render the methylation pattern as a string the length of the read.

    ``M`` marks a retained call, ``m`` one dropped for sitting in an excluded
    context, and ``.`` an unmodified base. By default every raw call is marked,
    including those below the probability threshold, which is how the original
    analysis was run. Set ``kept_only`` to mark only the calls that passed.
    """
    line = ["."] * length
    source = positions.kept if kept_only else positions.all_calls
    for position in source:
        if position < length:
            line[position] = "m" if position in positions.excluded else "M"
    if kept_only:
        for position in positions.excluded:
            if position < length:
                line[position] = "m"
    return "".join(line)
