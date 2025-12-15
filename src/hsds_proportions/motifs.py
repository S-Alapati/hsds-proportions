"""Motif definitions and the regexes built from them.

An hsdS allele is called from the recognition motif its methyltransferase
leaves behind. Each family is a pair of palindromic half sites separated by a
run of unconstrained bases, so the regex is assembled at run time once the
spacer length is known.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

#: Motif families of the WW2842 hsdS shufflon. Each entry gives the forward and
#: reverse strand recognition sequence with the spacer written as ``{spacer}``,
#: and the offsets within the match at which the adenine is expected to carry
#: the 6mA call.
DEFAULT_FAMILIES: dict[str, dict] = {
    "N1C1": {"fwd": "TCA{spacer}TGT", "rev": "ACA{spacer}TGA", "offsets": [2]},
    "N2C2": {"fwd": "GTA{spacer}TTA", "rev": "TAA{spacer}TAC", "offsets": [2]},
    "N1C2": {"fwd": "TCA{spacer}TTA", "rev": "TAA{spacer}TCA", "offsets": [2]},
    "N2C1": {"fwd": "GTA{spacer}TGT", "rev": "ACA{spacer}TAC", "offsets": [2]},
}

#: Sequence contexts whose 6mA calls are discarded before motifs are scored,
#: given as (context, offset of the adenine within the context). These are the
#: Dam and Dcm-like contexts that produce calls unrelated to the Type I system.
DEFAULT_EXCLUSIONS: list[tuple[str, int]] = [("CGCAG", 3), ("CCAGG", 2), ("GGACC", 2)]


@dataclass(frozen=True)
class MotifSet:
    """Compiled motif patterns for one spacer setting."""

    patterns: dict[str, tuple[re.Pattern, tuple[int, ...]]]
    exclusions: tuple[tuple[str, int], ...]
    min_spacer: int
    max_spacer: int
    families: tuple[str, ...] = field(default=())

    @property
    def family_names(self) -> tuple[str, ...]:
        return self.families


def build_motif_set(
    families: dict[str, dict] | None = None,
    min_spacer: int = 7,
    max_spacer: int = 7,
    exclusions: list[tuple[str, int]] | None = None,
) -> MotifSet:
    """Compile the forward and reverse pattern for every family.

    The spacer is a character class rather than a wildcard so that an N in the
    basecall never matches, which would otherwise inflate the counts.
    """
    families = families or DEFAULT_FAMILIES
    exclusions = DEFAULT_EXCLUSIONS if exclusions is None else exclusions

    if min_spacer < 0 or max_spacer < min_spacer:
        raise ValueError(f"invalid spacer range {min_spacer}-{max_spacer}")

    spacer = f"[ACGT]{{{min_spacer},{max_spacer}}}"
    patterns: dict[str, tuple[re.Pattern, tuple[int, ...]]] = {}

    for name, spec in families.items():
        missing = {"fwd", "rev"} - set(spec)
        if missing:
            raise ValueError(f"family {name} is missing {sorted(missing)}")
        offsets = tuple(spec.get("offsets", [2]))
        for strand in ("fwd", "rev"):
            label = f"{name}_{'Fwd' if strand == 'fwd' else 'Rev'}"
            patterns[label] = (re.compile(spec[strand].format(spacer=spacer)), offsets)

    return MotifSet(
        patterns=patterns,
        exclusions=tuple(tuple(e) for e in exclusions),
        min_spacer=min_spacer,
        max_spacer=max_spacer,
        families=tuple(families),
    )


def load_families(path: str | Path) -> dict[str, dict]:
    """Read motif families from JSON, so a different locus can be scored."""
    with open(path, encoding="utf-8") as handle:
        families = json.load(handle)
    if not isinstance(families, dict) or not families:
        raise ValueError(f"{path} does not contain a motif family mapping")
    return families


def parse_exclusion(text: str) -> tuple[str, int]:
    """Parse a ``SEQUENCE:OFFSET`` exclusion given on the command line."""
    sequence, _, offset = text.partition(":")
    sequence = sequence.strip().upper()
    if not sequence or not offset:
        raise ValueError(f"exclusion {text!r} should look like CCAGG:2")
    if set(sequence) - set("ACGT"):
        raise ValueError(f"exclusion {sequence!r} contains a non-ACGT base")
    index = int(offset)
    if not 0 <= index < len(sequence):
        raise ValueError(f"offset {index} falls outside {sequence!r}")
    return sequence, index
