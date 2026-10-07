"""Tests for modular TAD callers grammar, expansion, and consensus aggregation."""

import pytest
import sys
from pathlib import Path

# Add workflow/scripts to sys.path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "workflow" / "scripts"))

from callers import (
    expand_caller_configs,
    expand_caller_token,
    parse_caller_config,
    resolve_consensus_caller_configs,
)
from callers.consensustad import weighted_interval_scheduling


def test_expand_caller_token():
    assert expand_caller_token("topdom_w20") == ["topdom_w20"]
    assert expand_caller_token("topdom_w[20..23]") == [
        "topdom_w20",
        "topdom_w21",
        "topdom_w22",
        "topdom_w23",
    ]
    assert expand_caller_token("topdom_w[20..26:2]") == [
        "topdom_w20",
        "topdom_w22",
        "topdom_w24",
        "topdom_w26",
    ]


def test_expand_caller_configs():
    configs = ["topdom_w[9..11]", "spectraltad_l1", "cooltools_w250kb"]
    expanded = expand_caller_configs(configs)
    assert expanded == [
        "topdom_w9",
        "topdom_w10",
        "topdom_w11",
        "spectraltad_l1",
        "cooltools_w250kb",
    ]


def test_parse_caller_config_topdom():
    parsed = parse_caller_config("topdom_w20")
    assert parsed == {"caller": "topdom", "window_size": 20}


def test_parse_caller_config_spectraltad():
    parsed = parse_caller_config("spectraltad_l1")
    assert parsed == {"caller": "spectraltad", "levels": 1}


def test_parse_caller_config_cooltools():
    assert parse_caller_config("cooltools_w250kb") == {
        "caller": "cooltools",
        "window_bp": 250_000,
    }
    assert parse_caller_config("cooltools_w250000") == {
        "caller": "cooltools",
        "window_bp": 250_000,
    }
    assert parse_caller_config("cooltools_1mb") == {
        "caller": "cooltools",
        "window_bp": 1_000_000,
    }


def test_parse_caller_config_consensus_tokens():
    assert parse_caller_config("majority_vote") == {"caller": "majority_vote"}
    assert parse_caller_config("majority_vote_narrow") == {
        "caller": "majority_vote_narrow"
    }
    assert parse_caller_config("topdom_majority_vote") == {
        "caller": "topdom_majority_vote"
    }
    assert parse_caller_config("cooltools_majority_vote") == {
        "caller": "cooltools_majority_vote"
    }
    assert parse_caller_config("consensustad") == {"caller": "consensustad"}


def test_resolve_consensus_caller_configs():
    actual = ["topdom_w9", "topdom_w10", "topdom_w11", "cooltools_w3mb"]

    # "all" or "*"
    assert resolve_consensus_caller_configs("all", actual) == actual
    assert resolve_consensus_caller_configs("*", actual) == actual

    # Tool prefix selector
    assert resolve_consensus_caller_configs("topdom", actual) == [
        "topdom_w9",
        "topdom_w10",
        "topdom_w11",
    ]
    assert resolve_consensus_caller_configs("cooltools", actual) == ["cooltools_w3mb"]
    assert resolve_consensus_caller_configs("spectraltad", actual) == []

    # Explicit ranges / tokens
    assert resolve_consensus_caller_configs(["topdom_w[9..10]"], actual) == [
        "topdom_w9",
        "topdom_w10",
    ]
    assert resolve_consensus_caller_configs("topdom_w10", actual) == ["topdom_w10"]

    # Mixed list
    assert resolve_consensus_caller_configs(["topdom_w9", "cooltools"], actual) == [
        "topdom_w9",
        "cooltools_w3mb",
    ]

    # None or empty
    assert resolve_consensus_caller_configs(None, actual) == []
    assert resolve_consensus_caller_configs([], actual) == []


def test_parse_caller_config_invalid_prefix_error():
    # User's exact question: what happens if I write spectraltad_w1?
    with pytest.raises(ValueError, match="Caller 'spectraltad' requires parameter prefix 'l'"):
        parse_caller_config("spectraltad_w1")

    with pytest.raises(ValueError, match="Caller 'topdom' requires parameter prefix 'w'"):
        parse_caller_config("topdom_l1")

    with pytest.raises(ValueError, match="Unknown TAD caller 'unknown'"):
        parse_caller_config("unknown_w10")


def test_weighted_interval_scheduling_empty():
    assert weighted_interval_scheduling([]) == []


def test_weighted_interval_scheduling_disjoint():
    intervals = [
        {"tad_start": 0, "tad_stop": 100, "weight": 1.0},
        {"tad_start": 100, "tad_stop": 200, "weight": 1.0},
        {"tad_start": 200, "tad_stop": 300, "weight": 1.0},
    ]
    res = weighted_interval_scheduling(intervals)
    assert len(res) == 3
    assert [iv["tad_start"] for iv in res] == [0, 100, 200]


def test_weighted_interval_scheduling_overlapping():
    intervals = [
        {"tad_start": 0, "tad_stop": 150, "weight": 1.0},
        {"tad_start": 100, "tad_stop": 250, "weight": 3.0},  # Heavy interval overlapping both
        {"tad_start": 200, "tad_stop": 350, "weight": 1.0},
    ]
    res = weighted_interval_scheduling(intervals)
    assert len(res) == 1
    assert res[0]["tad_start"] == 100
    assert res[0]["weight"] == 3.0


def test_weighted_interval_scheduling_two_smaller_vs_one_larger():
    intervals = [
        {"tad_start": 0, "tad_stop": 100, "weight": 2.0},
        {"tad_start": 50, "tad_stop": 150, "weight": 2.5},  # Bridge interval
        {"tad_start": 100, "tad_stop": 200, "weight": 2.0},
    ]
    res = weighted_interval_scheduling(intervals)
    assert len(res) == 2
    assert [iv["tad_start"] for iv in res] == [0, 100]
