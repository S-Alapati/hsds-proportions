"""Motif definitions and the regexes built from them.

A Type I system recognises a bipartite site: two short half sites, each read by
one target recognition domain, separated by a run of unconstrained bases. Swap a
domain and the half site it reads changes, so the motif is what tells you which
allele a cell is carrying.

Motifs are written with IUPAC codes and ``{spacer}`` where the unconstrained run
falls, for example ``GTAY{spacer}TGT``. Nothing here is specific to any one
organism: the motifs come from a definition file you supply.
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import dataclass
from importlib import resources
from pathlib import Path

log = logging.getLogger(__name__)

#: IUPAC nucleotide codes and the bases each one stands for.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}

COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

SPACER_TOKEN = "{spacer}"


def reverse_complement(motif: str) -> str:
    """Reverse complement a motif, leaving the spacer token where it belongs."""
    if SPACER_TOKEN in motif:
        head, _, tail = motif.partition(SPACER_TOKEN)
        return f"{reverse_complement(tail)}{SPACER_TOKEN}{reverse_complement(head)}"
    return motif.upper().translate(COMPLEMENT)[::-1]


def to_regex(motif: str, min_spacer: int, max_spacer: int) -> str:
    """Turn an IUPAC motif into a regular expression.

    Every position becomes an explicit character class, so an N in the basecall
    never matches and cannot inflate the counts.
    """
    spacer = f"[ACGT]{{{min_spacer},{max_spacer}}}"
    parts = []
    # Split on the token before changing case, or the placeholder is mangled.
    for segment in motif.split(SPACER_TOKEN):
        rendered = ""
        for base in segment.upper():
            if base not in IUPAC:
                raise ValueError(f"{base!r} in {motif!r} is not an IUPAC code")
            bases = IUPAC[base]
            rendered += bases if len(bases) == 1 else f"[{bases}]"
        parts.append(rendered)
    return spacer.join(parts)


@dataclass(frozen=True)
class MotifSet:
    """Compiled patterns for one set of alleles."""

    patterns: dict[str, tuple[re.Pattern, tuple[int, ...]]]
    exclusions: tuple[tuple[str, int], ...]
    families: tuple[str, ...]
    spacers: dict[str, tuple[int, int]]
    name: str = "custom"

    @property
    def family_names(self) -> tuple[str, ...]:
        return self.families

    def describe(self) -> str:
        spans = {f"{lo}-{hi}" if lo != hi else str(lo) for lo, hi in self.spacers.values()}
        return f"{self.name}: {len(self.families)} alleles, spacer {', '.join(sorted(spans))}"


def build_motif_set(
    families: dict[str, dict],
    min_spacer: int = 7,
    max_spacer: int = 7,
    exclusions: list[tuple[str, int]] | None = None,
    name: str = "custom",
    check_complements: bool = True,
) -> MotifSet:
    """Compile the forward and reverse pattern for every allele.

    ``min_spacer`` and ``max_spacer`` apply to any family that does not set its
    own. A family may override them, which is needed when the two domains of a
    system space their half sites differently.
    """
    if min_spacer < 0 or max_spacer < min_spacer:
        raise ValueError(f"invalid spacer range {min_spacer}-{max_spacer}")

    patterns: dict[str, tuple[re.Pattern, tuple[int, ...]]] = {}
    spacers: dict[str, tuple[int, int]] = {}

    for family, spec in families.items():
        if "fwd" not in spec:
            raise ValueError(f"family {family} has no forward motif")
        forward = spec["fwd"]
        # A reverse motif can be given explicitly, or derived where the site is
        # read the same way on both strands, which is the usual case.
        reverse = spec.get("rev") or reverse_complement(forward)

        if check_complements and spec.get("rev"):
            expected = reverse_complement(forward)
            if expected.upper() != spec["rev"].upper():
                log.warning(
                    "family %s: reverse motif %s is not the reverse complement of %s "
                    "(expected %s). This is allowed, but check it is deliberate.",
                    family, spec["rev"], forward, expected,
                )

        low = int(spec.get("min_spacer", spec.get("spacer", min_spacer)))
        high = int(spec.get("max_spacer", spec.get("spacer", max_spacer)))
        if high < low:
            raise ValueError(f"family {family} has an invalid spacer range {low}-{high}")
        spacers[family] = (low, high)

        # A motif with no fixed base matches every read, which is almost always
        # an unedited skeleton rather than a real definition.
        fixed = [b for segment in forward.split(SPACER_TOKEN) for b in segment.upper() if b != "N"]
        if not fixed:
            raise ValueError(
                f"family {family} has no fixed base outside the spacer, so it would match "
                "every read. Replace the placeholder half sites with real motifs."
            )

        offsets = tuple(spec.get("offsets", [2]))
        for label, motif in ((f"{family}_Fwd", forward), (f"{family}_Rev", reverse)):
            patterns[label] = (re.compile(to_regex(motif, low, high)), offsets)

    return MotifSet(
        patterns=patterns,
        exclusions=tuple(tuple(e) for e in (exclusions or [])),
        families=tuple(families),
        spacers=spacers,
        name=name,
    )


# --- presets -------------------------------------------------------------

def _preset_dir():
    return resources.files(__package__) / "presets"


def list_presets() -> list[str]:
    """Names of the bundled presets."""
    return sorted(p.name[:-5] for p in _preset_dir().iterdir() if p.name.endswith(".json"))


def load_preset(name: str) -> dict:
    path = _preset_dir() / f"{name}.json"
    if not path.is_file():
        raise ValueError(f"unknown preset {name!r}, available: {', '.join(list_presets())}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_definition(path: str | Path) -> dict:
    """Read a motif definition from a JSON file."""
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict) or not data:
        raise ValueError(f"{path} does not contain a motif definition")
    # Allow either a bare mapping of families or a full definition document.
    return data if "families" in data else {"families": data}


def parse_exclusion(text: str) -> tuple[str, int]:
    """Parse a ``SEQUENCE:OFFSET`` exclusion given on the command line."""
    sequence, _, offset = text.partition(":")
    sequence = sequence.strip().upper()
    if not sequence or not offset:
        raise ValueError(f"exclusion {text!r} should look like CCAGG:2")
    unknown = set(sequence) - set(IUPAC)
    if unknown:
        raise ValueError(f"exclusion {sequence!r} contains {sorted(unknown)}")
    index = int(offset)
    if not 0 <= index < len(sequence):
        raise ValueError(f"offset {index} falls outside {sequence!r}")
    return sequence, index
