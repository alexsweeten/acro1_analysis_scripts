#!/usr/bin/env python3
"""Compare ACRO1 allele consensuses among six primates and with POLR2A CTD.

The primary peptide is a hypothetical translation of the oriented ACRO1 consensus
starting at nucleotide 1 (forward frame 0).  This phase was chosen because a
canonical ACRO1 unit is 147 nt (49 codons), while both named deletion alleles are
144 nt (48 codons).  All six frames are also evaluated as a sensitivity check.
"""

from __future__ import annotations

import gzip
import itertools
import math
import os
import random
import statistics
from collections import Counter
from pathlib import Path

import edlib
import pandas as pd


INPUT = Path(
    os.environ.get(
        "ACRO1_ALLELE_TABLE",
        Path(__file__).resolve().parents[1]
        / "inputs"
        / "all_primate_acro1_subfamily_locations.tsv.gz",
    )
)
OUTDIR = Path(
    os.environ.get("ACRO1_ALLELE_OUTPUT", Path(__file__).resolve().parents[1] / "results")
)

SPECIES_ORDER = [
    "hg002",
    "mPanTro3",
    "mPanPan1",
    "mGorGor1",
    "mPonAbe1",
    "mPonPyg2",
]
SPECIES_NAMES = {
    "hg002": "Human (HG002)",
    "mPanTro3": "Chimpanzee",
    "mPanPan1": "Bonobo",
    "mGorGor1": "Gorilla",
    "mPonAbe1": "Sumatran orangutan",
    "mPonPyg2": "Bornean orangutan",
}

ALLELE_ORDER = ["canonical_ACRO1", "internal_CCC", "boundary_terminal"]
ALLELE_NAMES = {
    "canonical_ACRO1": "Canonical ACRO1 (147 nt)",
    "internal_CCC": "Internal-CCC deletion (144 nt)",
    "boundary_terminal": "Terminal deletion (144 nt)",
}
EXPECTED_LENGTHS = {
    "canonical_ACRO1": 147,
    "internal_CCC": 144,
    "boundary_terminal": 144,
}

# Reviewed human POLR2A (UniProt P24928), residues 1728-1776 inclusive.
POLR2A_CTD_1728_1776 = "SPTSPSYSPTSPSYSPTSPNYSPTSPNYTPTSPSYSPTSPSYSPTSPNY"

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")


def translate(seq: str, frame: int = 0) -> str:
    seq = seq.upper()
    return "".join(
        CODON_TABLE.get(seq[i : i + 3], "X")
        for i in range(frame, len(seq) - 2, 3)
    )


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def consensus_and_support(sequences: list[str]) -> tuple[str, list[float], int]:
    lengths = {len(s) for s in sequences}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not equal length: {sorted(lengths)}")
    consensus = []
    support = []
    tied_positions = 0
    for column in zip(*sequences):
        counts = Counter(column)
        maximum = max(counts.values())
        winners = sorted(base for base, count in counts.items() if count == maximum)
        if len(winners) > 1:
            tied_positions += 1
        consensus.append(winners[0])
        support.append(maximum / len(column))
    return "".join(consensus), support, tied_positions


def species_balanced_consensus(
    species_consensuses: list[str], pooled_sequences: list[str]
) -> tuple[str, int, str]:
    """Plurality across species, with pooled counts only as deterministic tiebreaker."""
    consensus = []
    tied_positions = []
    for index, column in enumerate(zip(*species_consensuses), start=1):
        counts = Counter(column)
        maximum = max(counts.values())
        winners = sorted(base for base, count in counts.items() if count == maximum)
        if len(winners) > 1:
            tied_positions.append(index)
            pooled_counts = Counter(seq[index - 1] for seq in pooled_sequences)
            pooled_maximum = max(pooled_counts[base] for base in winners)
            winners = [base for base in winners if pooled_counts[base] == pooled_maximum]
        consensus.append(sorted(winners)[0])
    return (
        "".join(consensus),
        len(tied_positions),
        ",".join(str(position) for position in tied_positions),
    )


def global_alignment(a: str, b: str) -> dict[str, object]:
    result = edlib.align(a, b, mode="NW", task="path")
    nice = edlib.getNiceAlignment(result, a, b)
    aligned_a = nice["query_aligned"]
    aligned_b = nice["target_aligned"]
    matches = sum(x == y for x, y in zip(aligned_a, aligned_b))
    alignment_length = len(aligned_a)
    return {
        "matches": matches,
        "alignment_length": alignment_length,
        "identity": matches / alignment_length if alignment_length else math.nan,
        "edit_distance": int(result["editDistance"]),
        "aligned_a": aligned_a,
        "aligned_b": aligned_b,
        "match_line": nice["matched_aligned"],
    }


def composition(sequence: str) -> dict[str, float]:
    usable = [aa for aa in sequence if aa in AA_ORDER]
    denominator = len(usable)
    counts = Counter(usable)
    result = {aa: counts[aa] / denominator for aa in AA_ORDER}
    result.update(
        {
            "fraction_S": counts["S"] / denominator,
            "fraction_P": counts["P"] / denominator,
            "fraction_T": counts["T"] / denominator,
            "fraction_Y": counts["Y"] / denominator,
            "fraction_SPTY": sum(counts[aa] for aa in "SPTY") / denominator,
            "fraction_SPT": sum(counts[aa] for aa in "SPT") / denominator,
        }
    )
    return result


def composition_cosine(a: str, b: str) -> float:
    a_comp = composition(a)
    b_comp = composition(b)
    x = [a_comp[aa] for aa in AA_ORDER]
    y = [b_comp[aa] for aa in AA_ORDER]
    numerator = sum(i * j for i, j in zip(x, y))
    denominator = math.sqrt(sum(i * i for i in x) * sum(j * j for j in y))
    return numerator / denominator if denominator else math.nan


def peptide_metrics(peptide: str) -> dict[str, object]:
    comp = composition(peptide)
    return {
        "peptide_length_aa": len(peptide),
        "stop_count": peptide.count("*"),
        "S_count": peptide.count("S"),
        "P_count": peptide.count("P"),
        "T_count": peptide.count("T"),
        "Y_count": peptide.count("Y"),
        "SPTY_count": sum(peptide.count(aa) for aa in "SPTY"),
        "S_fraction": comp["fraction_S"],
        "P_fraction": comp["fraction_P"],
        "T_fraction": comp["fraction_T"],
        "Y_fraction": comp["fraction_Y"],
        "SPTY_fraction": comp["fraction_SPTY"],
    }


def ctd_metrics(peptide: str) -> dict[str, object]:
    alignment = global_alignment(peptide, POLR2A_CTD_1728_1776)
    longest, peptide_start, ctd_start = longest_common_substring(
        peptide, POLR2A_CTD_1728_1776
    )
    return {
        "ctd_matches": alignment["matches"],
        "ctd_alignment_length": alignment["alignment_length"],
        "ctd_global_identity": alignment["identity"],
        "ctd_edit_distance": alignment["edit_distance"],
        "ctd_composition_cosine": composition_cosine(peptide, POLR2A_CTD_1728_1776),
        "ctd_longest_exact_substring": longest,
        "ctd_longest_exact_substring_length": len(longest),
        "ctd_longest_exact_substring_peptide_start_1based": peptide_start,
        "ctd_longest_exact_substring_ctd_start_1based": ctd_start,
        "ctd_aligned_peptide": alignment["aligned_a"],
        "ctd_match_line": alignment["match_line"],
        "ctd_aligned_reference": alignment["aligned_b"],
    }


def longest_common_substring(a: str, b: str) -> tuple[str, int, int]:
    best = ""
    best_a = 0
    best_b = 0
    for index_a in range(len(a)):
        for index_b in range(len(b)):
            length = 0
            while (
                index_a + length < len(a)
                and index_b + length < len(b)
                and a[index_a + length] == b[index_b + length]
            ):
                length += 1
            if length > len(best):
                best = a[index_a : index_a + length]
                best_a = index_a + 1
                best_b = index_b + 1
    return best, best_a, best_b


def alignment_identity_only(a: str, b: str) -> float:
    result = edlib.align(a, b, mode="NW", task="path")
    nice = edlib.getNiceAlignment(result, a, b)
    aligned_a = nice["query_aligned"]
    aligned_b = nice["target_aligned"]
    return sum(x == y for x, y in zip(aligned_a, aligned_b)) / len(aligned_a)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    with gzip.open(INPUT, "rt") as handle:
        calls = pd.read_csv(handle, sep="\t", low_memory=False)
    calls = calls[
        calls["species_id"].isin(SPECIES_ORDER)
        & calls["subfamily_label"].isin(ALLELE_ORDER)
    ].copy()
    calls["sequence"] = calls["sequence"].str.upper()

    invalid = calls[~calls["sequence"].str.fullmatch("[ACGT]+")]
    if not invalid.empty:
        raise ValueError(f"Found {len(invalid)} sequences with non-ACGT characters")
    for allele, expected_length in EXPECTED_LENGTHS.items():
        observed = set(calls.loc[calls["subfamily_label"] == allele, "sequence"].str.len())
        if observed != {expected_length}:
            raise ValueError(f"{allele}: expected {expected_length}, observed {observed}")

    per_species_rows = []
    nucleotide_fasta = []
    peptide_fasta = []
    species_consensus_lookup: dict[tuple[str, str], str] = {}
    for species_id in SPECIES_ORDER:
        for allele in ALLELE_ORDER:
            subset = calls[
                (calls["species_id"] == species_id)
                & (calls["subfamily_label"] == allele)
            ]
            sequences = subset["sequence"].tolist()
            if not sequences:
                raise ValueError(f"No sequences for {species_id} / {allele}")
            consensus, support, ties = consensus_and_support(sequences)
            peptide = translate(consensus)
            species_consensus_lookup[(species_id, allele)] = consensus
            row = {
                "species_id": species_id,
                "species_name": SPECIES_NAMES[species_id],
                "allele": allele,
                "allele_name": ALLELE_NAMES[allele],
                "unit_count": len(sequences),
                "nucleotide_length": len(consensus),
                "consensus_nt": consensus,
                "display_nt_bases_1_25_and_125_end": f"{consensus[:25]}...{consensus[124:]}",
                "mean_modal_base_fraction": sum(support) / len(support),
                "minimum_modal_base_fraction": min(support),
                "positions_modal_fraction_ge_0.5": sum(x >= 0.5 for x in support),
                "tied_consensus_positions": ties,
                "forward_frame0_peptide": peptide,
                **peptide_metrics(peptide),
                **ctd_metrics(peptide),
            }
            per_species_rows.append(row)
            nucleotide_fasta.append((f"{species_id}|{allele}|n={len(sequences)}", consensus))
            peptide_fasta.append((f"{species_id}|{allele}|frame=+0|n={len(sequences)}", peptide))

    per_species = pd.DataFrame(per_species_rows)
    per_species.to_csv(OUTDIR / "species_allele_consensus_and_ctd.tsv", sep="\t", index=False)
    write_fasta(OUTDIR / "species_allele_consensus_nt.fasta", nucleotide_fasta)
    write_fasta(OUTDIR / "species_allele_frame0_peptides.fasta", peptide_fasta)

    pairwise_rows = []
    for allele in ALLELE_ORDER:
        for species_a, species_b in itertools.combinations(SPECIES_ORDER, 2):
            nt_a = species_consensus_lookup[(species_a, allele)]
            nt_b = species_consensus_lookup[(species_b, allele)]
            aa_a = translate(nt_a)
            aa_b = translate(nt_b)
            nt_alignment = global_alignment(nt_a, nt_b)
            aa_alignment = global_alignment(aa_a, aa_b)
            pairwise_rows.append(
                {
                    "allele": allele,
                    "species_a": species_a,
                    "species_b": species_b,
                    "nt_matches": nt_alignment["matches"],
                    "nt_alignment_length": nt_alignment["alignment_length"],
                    "nt_identity": nt_alignment["identity"],
                    "aa_matches": aa_alignment["matches"],
                    "aa_alignment_length": aa_alignment["alignment_length"],
                    "aa_identity": aa_alignment["identity"],
                }
            )
    pairwise = pd.DataFrame(pairwise_rows)
    pairwise.to_csv(OUTDIR / "within_allele_pairwise_primate_identity.tsv", sep="\t", index=False)

    conservation_rows = []
    for allele in ALLELE_ORDER:
        allele_pairwise = pairwise[pairwise["allele"] == allele]
        consensuses = [species_consensus_lookup[(species, allele)] for species in SPECIES_ORDER]
        peptides = [translate(sequence) for sequence in consensuses]
        variable_nt_positions = sum(len(set(column)) > 1 for column in zip(*consensuses))
        variable_aa_positions = sum(len(set(column)) > 1 for column in zip(*peptides))
        conservation_rows.append(
            {
                "allele": allele,
                "allele_name": ALLELE_NAMES[allele],
                "all_six_nt_consensuses_identical": len(set(consensuses)) == 1,
                "all_six_frame0_peptides_identical": len(set(peptides)) == 1,
                "variable_nt_positions_among_six_consensuses": variable_nt_positions,
                "variable_aa_positions_among_six_consensuses": variable_aa_positions,
                "pairwise_nt_identity_min": allele_pairwise["nt_identity"].min(),
                "pairwise_nt_identity_mean": allele_pairwise["nt_identity"].mean(),
                "pairwise_nt_identity_max": allele_pairwise["nt_identity"].max(),
                "pairwise_aa_identity_min": allele_pairwise["aa_identity"].min(),
                "pairwise_aa_identity_mean": allele_pairwise["aa_identity"].mean(),
                "pairwise_aa_identity_max": allele_pairwise["aa_identity"].max(),
                "unique_nt_consensus_count": len(set(consensuses)),
                "unique_frame0_peptide_count": len(set(peptides)),
            }
        )
    conservation = pd.DataFrame(conservation_rows)
    conservation.to_csv(OUTDIR / "allele_cross_primate_conservation_summary.tsv", sep="\t", index=False)

    variable_rows = []
    for allele in ALLELE_ORDER:
        allele_consensuses = {
            species: species_consensus_lookup[(species, allele)] for species in SPECIES_ORDER
        }
        allele_peptides = {
            species: translate(allele_consensuses[species]) for species in SPECIES_ORDER
        }
        for molecule, sequences in [("nucleotide", allele_consensuses), ("amino_acid", allele_peptides)]:
            for position, column in enumerate(zip(*(sequences[s] for s in SPECIES_ORDER)), start=1):
                if len(set(column)) == 1:
                    continue
                row = {
                    "allele": allele,
                    "molecule": molecule,
                    "position_1based_within_native_allele": position,
                    "distinct_states": ",".join(sorted(set(column))),
                }
                row.update({species: sequences[species][position - 1] for species in SPECIES_ORDER})
                variable_rows.append(row)
    pd.DataFrame(variable_rows).to_csv(
        OUTDIR / "cross_primate_variable_consensus_positions.tsv", sep="\t", index=False
    )

    balanced_rows = []
    balanced_nt_fasta = []
    balanced_aa_fasta = []
    alignment_blocks = []
    rng = random.Random(17281776)
    shuffled_ctd_sequences = []
    for _ in range(10_000):
        residues = list(POLR2A_CTD_1728_1776)
        rng.shuffle(residues)
        shuffled_ctd_sequences.append("".join(residues))
    for allele in ALLELE_ORDER:
        per_species_consensuses = [
            species_consensus_lookup[(species, allele)] for species in SPECIES_ORDER
        ]
        pooled_sequences = calls.loc[calls["subfamily_label"] == allele, "sequence"].tolist()
        consensus, tie_count, tie_positions = species_balanced_consensus(
            per_species_consensuses, pooled_sequences
        )
        peptide = translate(consensus)
        ctd = ctd_metrics(peptide)
        shuffled_identities = [
            alignment_identity_only(peptide, shuffled_ctd)
            for shuffled_ctd in shuffled_ctd_sequences
        ]
        observed_identity = float(ctd["ctd_global_identity"])
        exceedances = sum(identity >= observed_identity for identity in shuffled_identities)
        balanced_rows.append(
            {
                "allele": allele,
                "allele_name": ALLELE_NAMES[allele],
                "species_count": len(SPECIES_ORDER),
                "total_units": len(pooled_sequences),
                "nucleotide_length": len(consensus),
                "species_balanced_consensus_nt": consensus,
                "species_plurality_tied_position_count": tie_count,
                "species_plurality_tied_positions_1based": tie_positions,
                "forward_frame0_peptide": peptide,
                **peptide_metrics(peptide),
                **ctd,
                "ctd_order_shuffle_count": len(shuffled_identities),
                "ctd_order_shuffle_mean_identity": statistics.mean(shuffled_identities),
                "ctd_order_shuffle_sd_identity": statistics.pstdev(shuffled_identities),
                "ctd_order_shuffle_empirical_p_ge_observed_plus1": (
                    exceedances + 1
                )
                / (len(shuffled_identities) + 1),
            }
        )
        balanced_nt_fasta.append((f"six_primates_species_balanced|{allele}", consensus))
        balanced_aa_fasta.append((f"six_primates_species_balanced|{allele}|frame=+0", peptide))
        alignment_blocks.extend(
            [
                ALLELE_NAMES[allele],
                f"ACRO1  {ctd['ctd_aligned_peptide']}",
                f"        {ctd['ctd_match_line']}",
                f"POLR2A {ctd['ctd_aligned_reference']}",
                "",
            ]
        )
    balanced = pd.DataFrame(balanced_rows)
    balanced.to_csv(OUTDIR / "species_balanced_allele_consensus_and_ctd.tsv", sep="\t", index=False)
    write_fasta(OUTDIR / "species_balanced_allele_consensus_nt.fasta", balanced_nt_fasta)
    write_fasta(OUTDIR / "species_balanced_allele_frame0_peptides.fasta", balanced_aa_fasta)
    (OUTDIR / "species_balanced_allele_vs_polr2a_alignments.txt").write_text(
        "POLR2A P24928 residues 1728-1776 (1-based, inclusive)\n"
        f"Reference: {POLR2A_CTD_1728_1776}\n\n" + "\n".join(alignment_blocks)
    )

    balanced_pairwise_rows = []
    for allele_a, allele_b in itertools.combinations(ALLELE_ORDER, 2):
        a = balanced.loc[balanced["allele"] == allele_a].iloc[0]
        b = balanced.loc[balanced["allele"] == allele_b].iloc[0]
        for molecule, column in [
            ("nucleotide", "species_balanced_consensus_nt"),
            ("amino_acid", "forward_frame0_peptide"),
        ]:
            alignment = global_alignment(a[column], b[column])
            balanced_pairwise_rows.append(
                {
                    "allele_a": allele_a,
                    "allele_b": allele_b,
                    "molecule": molecule,
                    "matches": alignment["matches"],
                    "alignment_length": alignment["alignment_length"],
                    "global_identity": alignment["identity"],
                    "edit_distance": alignment["edit_distance"],
                    "aligned_a": alignment["aligned_a"],
                    "match_line": alignment["match_line"],
                    "aligned_b": alignment["aligned_b"],
                }
            )
    pd.DataFrame(balanced_pairwise_rows).to_csv(
        OUTDIR / "species_balanced_allele_pairwise_identity.tsv", sep="\t", index=False
    )

    frame_rows = []
    for row in per_species_rows:
        nt = row["consensus_nt"]
        for strand, strand_sequence in [("+", nt), ("-", reverse_complement(nt))]:
            for frame in range(3):
                peptide = translate(strand_sequence, frame)
                ctd = ctd_metrics(peptide)
                frame_rows.append(
                    {
                        "species_id": row["species_id"],
                        "species_name": row["species_name"],
                        "allele": row["allele"],
                        "strand": strand,
                        "frame": frame,
                        "frame_label": f"{'forward' if strand == '+' else 'reverse'}_{frame}",
                        "peptide": peptide,
                        **peptide_metrics(peptide),
                        **ctd,
                    }
                )
    frames = pd.DataFrame(frame_rows)
    frames["best_identity_frame_for_species_allele"] = False
    best_indices = frames.groupby(["species_id", "allele"])["ctd_global_identity"].idxmax()
    frames.loc[best_indices, "best_identity_frame_for_species_allele"] = True
    frames.to_csv(OUTDIR / "all_six_frames_vs_polr2a_ctd.tsv", sep="\t", index=False)

    ctd_composition = peptide_metrics(POLR2A_CTD_1728_1776)
    ctd_summary = pd.DataFrame(
        [
            {
                "reference": "UniProt P24928 human POLR2A",
                "residues_1based_inclusive": "1728-1776",
                "sequence": POLR2A_CTD_1728_1776,
                **ctd_composition,
            }
        ]
    )
    ctd_summary.to_csv(OUTDIR / "polr2a_ctd_reference.tsv", sep="\t", index=False)

    print(f"Read {len(calls):,} strict ACRO1 calls")
    print(f"Wrote results to {OUTDIR}")
    print("\nCross-primate conservation:")
    print(conservation.to_string(index=False))
    print("\nSpecies-balanced allele products vs POLR2A CTD:")
    print(
        balanced[
            [
                "allele",
                "forward_frame0_peptide",
                "stop_count",
                "S_fraction",
                "P_fraction",
                "T_fraction",
                "SPTY_fraction",
                "ctd_matches",
                "ctd_alignment_length",
                "ctd_global_identity",
                "ctd_composition_cosine",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
