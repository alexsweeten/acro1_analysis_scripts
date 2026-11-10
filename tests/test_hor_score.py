#!/usr/bin/env python3
"""Unit tests for ACRE-scale HOR partner counting."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "analyses"
    / "12_hor_score"
    / "scripts"
    / "compute_hor_score.py"
)
SPEC = importlib.util.spec_from_file_location("compute_hor_score", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
HOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HOR)


class ScoreDiagonalRunsTest(unittest.TestCase):
    def test_distinct_bidirectional_partner_counting(self) -> None:
        # Offset 1 has a four-pair run followed by a failure; offset 2 has a
        # four-pair run. Partners repeated across blocks would still be counted
        # only once because the implementation uses sets.
        variant_offsets = [
            [0, 0, 0, 0, 100],
            [0, 0, 0, 0],
            [100, 100, 100],
            [100, 100],
            [100],
        ]
        clusters = ["H1", "H1", "H1", "H1", "H2", "H2"]
        partners, same_cluster, block_lengths = HOR.score_diagonal_runs(
            variant_offsets,
            mismatch_threshold=0,
            clusters=clusters,
            min_block_length=4,
        )
        self.assertEqual(partners, [2, 3, 4, 4, 2, 1])
        self.assertEqual(same_cluster, [2, 3, 3, 2, 0, 0])
        self.assertEqual(block_lengths, [4, 4])

    def test_minimum_block_length_is_enforced(self) -> None:
        partners, same_cluster, block_lengths = HOR.score_diagonal_runs(
            [[0, 0, 0, 0]],
            mismatch_threshold=0,
            clusters=["H1"] * 5,
            min_block_length=5,
        )
        self.assertEqual(partners, [0, 0, 0, 0, 0])
        self.assertEqual(same_cluster, [0, 0, 0, 0, 0])
        self.assertEqual(block_lengths, [])


if __name__ == "__main__":
    unittest.main()
