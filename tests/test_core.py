"""Tests for the parts that do not need pysam or external tools."""

from __future__ import annotations

import pytest

from hsds_proportions.assign import AMBIGUOUS, UNASSIGNED, MotifHit, assign_allele, find_motifs
from hsds_proportions.modifications import annotation_string, filter_calls, parse_mm_ml
from hsds_proportions.motifs import build_motif_set, parse_exclusion
from hsds_proportions.summary import ErrorModel, proportions


def test_spacer_length_is_respected():
    motifs = build_motif_set(min_spacer=7, max_spacer=7)
    pattern, _ = motifs.patterns["N1C1_Fwd"]
    assert pattern.search("TCA" + "A" * 7 + "TGT")
    assert not pattern.search("TCA" + "A" * 6 + "TGT")


def test_spacer_will_not_match_an_ambiguous_base():
    motifs = build_motif_set()
    pattern, _ = motifs.patterns["N1C1_Fwd"]
    assert not pattern.search("TCA" + "A" * 6 + "N" + "TGT")


def test_motif_is_only_reported_when_the_adenine_is_methylated():
    motifs = build_motif_set()
    sequence = "GG" + "TCA" + "ACGTACG" + "TGT" + "GG"
    annotation = "." * len(sequence)
    assert find_motifs(sequence, set(), annotation, motifs) == []
    hits = find_motifs(sequence, {4}, annotation, motifs)
    assert [h.family for h in hits] == ["N1C1"]
    assert hits[0].sequence == "TCAACGTACGTGT"


def test_a_read_with_one_family_is_assigned_to_it():
    families = ("N1C1", "N2C2", "N1C2", "N2C1")
    hits = [MotifHit("N1C1", "N1C1_Fwd", "", "")] * 3
    assert assign_allele(hits, families) == "N1C1"


def test_a_clear_majority_wins_once_there_are_enough_motifs():
    families = ("N1C1", "N2C2", "N1C2", "N2C1")
    hits = [MotifHit("N1C1", "N1C1_Fwd", "", "")] * 9 + [MotifHit("N2C2", "N2C2_Fwd", "", "")]
    assert assign_allele(hits, families) == "N1C1"


def test_a_split_read_is_held_as_ambiguous():
    families = ("N1C1", "N2C2", "N1C2", "N2C1")
    hits = [MotifHit("N1C1", "N1C1_Fwd", "", "")] * 3 + [MotifHit("N2C2", "N2C2_Fwd", "", "")] * 3
    assert assign_allele(hits, families) == AMBIGUOUS


def test_a_read_without_motifs_is_unassigned():
    assert assign_allele([], ("N1C1",)) == UNASSIGNED


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
    shares = proportions({"N1C1": 80, "N2C2": 20})
    assert shares == {"N1C1": 80.0, "N2C2": 20.0}


def test_proportions_of_nothing_are_zero():
    assert proportions({"N1C1": 0, "N2C2": 0}) == {"N1C1": 0.0, "N2C2": 0.0}


def test_error_model_reproduces_the_published_value():
    # Q99 with a probability floor of 255 gave 4.9961% in the WW2842 analysis.
    rate = ErrorModel().error_rate(mean_accuracy_pct=99.0, min_probability=255.0)
    assert rate == pytest.approx(4.9961, abs=1e-4)


def test_exclusions_are_parsed_and_validated():
    assert parse_exclusion("CCAGG:2") == ("CCAGG", 2)
    with pytest.raises(ValueError):
        parse_exclusion("CCAGG")
    with pytest.raises(ValueError):
        parse_exclusion("CCAGG:9")
    with pytest.raises(ValueError):
        parse_exclusion("CCAXG:2")
