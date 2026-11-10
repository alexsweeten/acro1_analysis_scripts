#!/usr/bin/env python3
"""Conservative species/assembly-level null test for ACRO1 consensus ORFs."""

from __future__ import annotations

import csv
import random
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "primate_acro1_orf_test"
STOP = {"TAA", "TAG", "TGA"}
RC = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(RC)[::-1]


def metrics(sequence: str) -> tuple[int, int, int]:
    longest_stopfree = 0
    longest_atg = 0
    zero_frames = 0
    for oriented in (sequence, reverse_complement(sequence)):
        for frame in range(3):
            codons = [oriented[index:index + 3] for index in range(frame, len(oriented) - 2, 3)]
            stops = [index for index, codon in enumerate(codons) if codon in STOP]
            if not stops:
                zero_frames += 1
            boundaries = [-1] + stops + [len(codons)]
            longest_stopfree = max(longest_stopfree, max(right - left - 1 for left, right in zip(boundaries, boundaries[1:])))
            active = []
            for index, codon in enumerate(codons):
                if codon == "ATG":
                    active.append(index)
                if codon in STOP:
                    if active:
                        longest_atg = max(longest_atg, index - min(active))
                    active = []
            if active:
                longest_atg = max(longest_atg, len(codons) - min(active))
    return longest_atg, longest_stopfree, zero_frames


def dinucleotide_shuffle(sequence: str, rng: random.Random) -> str:
    adjacency: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(sequence, sequence[1:]):
        adjacency[left].append(right)
    for edges in adjacency.values():
        rng.shuffle(edges)
    stack = [sequence[0]]
    trail = []
    while stack:
        node = stack[-1]
        if adjacency[node]:
            stack.append(adjacency[node].pop())
        else:
            trail.append(stack.pop())
    return "".join(reversed(trail))


def empirical_p(observed: int, null: list[int]) -> float:
    return (1 + sum(value >= observed for value in null)) / (1 + len(null))


permutations = 999
rows = list(csv.DictReader((OUT / "assembly_monomer_consensuses.tsv").open(), delimiter="\t"))
results = []
for row_index, row in enumerate(rows):
    sequence = row["consensus_sequence"]
    observed_atg, observed_stopfree, observed_zero = metrics(sequence * 10)
    candidate_frame_stops = sum(
        sequence[index:index + 3] in STOP for index in range(0, len(sequence) - 2, 3)
    )
    rng = random.Random(500 + row_index)
    null_atg = []
    null_stopfree = []
    for permutation in range(permutations):
        shuffled = dinucleotide_shuffle(sequence, rng)
        atg, stopfree, zero = metrics(shuffled * 10)
        null_atg.append(atg)
        null_stopfree.append(stopfree)
    results.append({
        "key": row["key"],
        "species": row["species"],
        "consensus_length": len(sequence),
        "candidate_plus0_stops_single_monomer": candidate_frame_stops,
        "observed_max_atg_open_aa_10x": observed_atg,
        "null_mean_max_atg_open_aa_10x": sum(null_atg) / permutations,
        "p_max_atg_open_aa": empirical_p(observed_atg, null_atg),
        "observed_max_stopfree_aa_10x": observed_stopfree,
        "null_mean_max_stopfree_aa_10x": sum(null_stopfree) / permutations,
        "p_max_stopfree_aa": empirical_p(observed_stopfree, null_stopfree),
        "observed_zero_stop_frames_10x": observed_zero,
    })

with (OUT / "assembly_consensus_null_test.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(results[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(results)

print(OUT / "assembly_consensus_null_test.tsv")
