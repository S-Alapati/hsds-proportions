"""Tests for the parts that do not need pysam or external tools."""

from __future__ import annotations

import pytest

from hsds_proportions.assign import AMBIGUOUS, UNASSIGNED, MotifHit, assign_allele, find_motifs
from hsds_proportions.modifications import annotation_string, filter_calls, parse_mm_ml
from pathlib import Path

from hsds_proportions.motifs import (
    build_motif_set,
    list_presets,
    load_definition,
    load_preset,
    parse_exclusion,
    reverse_complement,
    to_regex,
)
from hsds_proportions.summary import ErrorModel, proportions


EXAMPLE = {
    "alleleA": {"fwd": "TCA{spacer}TGT", "rev": "ACA{spacer}TGA", "spacer": 7},
    "alleleB": {"fwd": "GTAY{spacer}TTA", "rev": "TAA{spacer}RTAC", "spacer": 6},
}


def test_spacer_length_is_respected():
    motifs = build_motif_set(EXAMPLE, min_spacer=7, max_spacer=7)
    pattern, _ = motifs.patterns["alleleA_Fwd"]
    assert pattern.search("TCA" + "A" * 7 + "TGT")
    assert not pattern.search("TCA" + "A" * 6 + "TGT")


def test_spacer_will_not_match_an_ambiguous_base():
    motifs = build_motif_set(EXAMPLE)
    pattern, _ = motifs.patterns["alleleA_Fwd"]
    assert not pattern.search("TCA" + "A" * 6 + "N" + "TGT")


def test_motif_is_only_reported_when_the_adenine_is_methylated():
    motifs = build_motif_set(EXAMPLE)
    sequence = "GG" + "TCA" + "ACGTACG" + "TGT" + "GG"
    annotation = "." * len(sequence)
    assert find_motifs(sequence, set(), annotation, motifs) == []
    hits = find_motifs(sequence, {4}, annotation, motifs)
    assert [h.family for h in hits] == ["alleleA"]
    assert hits[0].sequence == "TCAACGTACGTGT"


FAMILIES = ("alleleA", "alleleB", "alleleC", "alleleD")


def _hits(family, count):
    return [MotifHit(family, f"{family}_Fwd", "", "")] * count


def test_a_read_with_one_family_is_assigned_to_it():
    assert assign_allele(_hits("alleleA", 3), FAMILIES) == "alleleA"


def test_a_clear_majority_wins_once_there_are_enough_motifs():
    hits = _hits("alleleA", 9) + _hits("alleleB", 1)
    assert assign_allele(hits, FAMILIES) == "alleleA"


def test_a_split_read_is_held_as_ambiguous():
    hits = _hits("alleleA", 3) + _hits("alleleB", 3)
    assert assign_allele(hits, FAMILIES) == AMBIGUOUS


def test_a_read_without_motifs_is_unassigned():
    assert assign_allele([], FAMILIES) == UNASSIGNED


def test_mm_deltas_land_on_the_right_adenines():
    # Adenines sit at 0, 3 and 6. A delta of 0 takes the first, then 1 skips
    # the second and takes the third.
    sequence = "ACCACCACC"
    calls = parse_mm_ml(sequence, "A+a.,0,1;", [200, 100])
    assert calls == {0: ("A+a.", 200), 6: ("A+a.", 100)}


def test_a_short_ml_tag_does_not_raise():
    assert parse_mm_ml("ACCACCACC", "A+a.,0,1;", [200]) == {}


def test_calls_in_an_excluded_context_are_dropped():
    sequence = "TTCCAGGTT"          # CCAGG with the adenine at offset 2, position 4
    calls = {4: ("A+a.", 255)}
    result = filter_calls(sequence, calls, "A+a.", 255, (("CCAGG", 2),))
    assert result.kept == set()
    assert result.excluded == {4}
    assert result.all_calls == {4}


def test_calls_below_the_threshold_are_not_kept():
    sequence = "TTTATTT"
    result = filter_calls(sequence, {3: ("A+a.", 200)}, "A+a.", 255, ())
    assert result.kept == set()
    assert result.all_calls == {3}


def test_annotation_marks_kept_and_excluded_differently():
    sequence = "TTCCAGGTA"
    calls = {4: ("A+a.", 255), 8: ("A+a.", 255)}
    positions = filter_calls(sequence, calls, "A+a.", 255, (("CCAGG", 2),))
    assert annotation_string(len(sequence), positions) == "....m...M"
    assert annotation_string(len(sequence), positions, kept_only=True) == "....m...M"


def test_proportions_sum_to_one_hundred():
    shares = proportions({"alleleA": 80, "alleleB": 20})
    assert shares == {"alleleA": 80.0, "alleleB": 20.0}


def test_proportions_of_nothing_are_zero():
    assert proportions({"alleleA": 0, "alleleB": 0}) == {"alleleA": 0.0, "alleleB": 0.0}


def test_the_error_model_multiplies_its_three_terms():
    # 1 - (0.999 * 0.99^5 * 1.0), as a percentage.
    rate = ErrorModel().error_rate(mean_accuracy_pct=99.0, min_probability=255.0)
    assert rate == pytest.approx(4.9961, abs=1e-4)


def test_a_lower_probability_floor_raises_the_error_rate():
    strict = ErrorModel().error_rate(99.0, 255.0)
    loose = ErrorModel().error_rate(99.0, 200.0)
    assert loose > strict


def test_exclusions_are_parsed_and_validated():
    assert parse_exclusion("CCAGG:2") == ("CCAGG", 2)
    with pytest.raises(ValueError):
        parse_exclusion("CCAGG")
    with pytest.raises(ValueError):
        parse_exclusion("CCAGG:9")
    with pytest.raises(ValueError):
        parse_exclusion("CCAXG:2")


def test_iupac_codes_become_character_classes():
    assert to_regex("GTAY{spacer}TGT", 6, 6) == "GTA[CT][ACGT]{6,6}TGT"
    assert to_regex("ACA{spacer}RTAC", 6, 6) == "ACA[ACGT]{6,6}[AG]TAC"


def test_a_degenerate_position_is_enforced():
    motifs = build_motif_set({"alleleX": {"fwd": "GTAY{spacer}TGT", "spacer": 6}})
    pattern, _ = motifs.patterns["alleleX_Fwd"]
    assert pattern.search("GTAC" + "A" * 6 + "TGT")
    assert pattern.search("GTAT" + "A" * 6 + "TGT")
    assert not pattern.search("GTAA" + "A" * 6 + "TGT")


def test_reverse_complement_keeps_the_spacer_in_place():
    assert reverse_complement("GTAY{spacer}TGT") == "ACA{spacer}RTAC"
    assert reverse_complement("TCA{spacer}TTA") == "TAA{spacer}TGA"


def test_a_missing_reverse_motif_is_derived():
    motifs = build_motif_set({"X": {"fwd": "TCA{spacer}TTA", "spacer": 7}})
    pattern, _ = motifs.patterns["X_Rev"]
    assert pattern.search("TAA" + "C" * 7 + "TGA")


def test_families_may_set_their_own_spacer():
    motifs = build_motif_set(EXAMPLE)
    assert motifs.spacers == {"alleleA": (7, 7), "alleleB": (6, 6)}


def test_the_worked_example_compiles():
    definition = load_definition(Path(__file__).parent.parent / "examples" / "four_allele_example.json")
    motif_set = build_motif_set(
        definition["families"],
        exclusions=[tuple(e) for e in definition.get("exclusions", [])],
        name="example",
    )
    assert len(motif_set.family_names) == 4
    assert motif_set.spacers["alleleC"] == (6, 6)


def test_a_motif_of_only_ns_is_rejected():
    with pytest.raises(ValueError, match="no fixed base"):
        build_motif_set({"skeleton": {"fwd": "NNN{spacer}NNN", "spacer": 6}})


def test_the_bundled_template_is_a_skeleton_and_will_not_run():
    with pytest.raises(ValueError, match="no fixed base"):
        build_motif_set(load_preset("template")["families"])


def test_a_preset_is_bundled_for_copying():
    assert "template" in list_presets()
