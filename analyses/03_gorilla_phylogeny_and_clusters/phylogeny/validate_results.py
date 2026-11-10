#!/usr/bin/env python3
"""Validate all-gorilla extraction, clustering, alignment, tree, and figures."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
EXPECTED_TIPS = 1232


def fasta_dimensions(path: Path) -> tuple[list[str], int | None]:
    names: list[str] = []
    lengths: list[int] = []
    current = 0
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if names:
                    lengths.append(current)
                names.append(line[1:].strip().split()[0])
                current = 0
            else:
                current += len(line.strip())
    if names:
        lengths.append(current)
    if len(set(lengths)) > 1:
        return names, None
    return names, lengths[0] if lengths else None


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def main() -> None:
    metadata = pd.read_csv(HERE / "gorilla_tip_metadata.tsv", sep="\t")
    base = pd.read_csv(HERE / "gorilla_tip_metadata.base.tsv", sep="\t")
    clusters = pd.read_csv(HERE / "gorilla_de_novo_clusters.tsv", sep="\t")
    parameters = json.loads((HERE / "cluster_sweep_parameters.json").read_text())
    expected = set(metadata["tip_id"])
    assert len(metadata) == EXPECTED_TIPS
    assert metadata["tip_id"].nunique() == EXPECTED_TIPS
    assert len(clusters) == EXPECTED_TIPS
    assert not clusters.isna().any().any()
    assert clusters["de_novo_cluster"].nunique() == parameters["best_k"] == 6
    assert metadata["synteny_group"].nunique() == 11
    assert not {"cluster8", "best_k", "best_cluster", "cluster8_centroid_distance"} & set(base.columns)

    fasta_files = [
        "gorilla_acres.fasta",
        "gorilla_acres.aligned.raw.fasta",
        "gorilla_acres.aligned.fasta",
        "gorilla_acres.aligned.occupancy02.fasta",
        "gorilla_acres.aligned.occupancy05.fasta",
        "gorilla_acres.aligned.occupancy10.fasta",
    ]
    for filename in fasta_files:
        names, width = fasta_dimensions(HERE / filename)
        normalized = {name[3:] if name.startswith("_R_") else name for name in names}
        assert len(names) == EXPECTED_TIPS, (filename, len(names))
        assert len(normalized) == EXPECTED_TIPS
        assert normalized == expected
        if ".aligned" in filename:
            assert width is not None and width > 0, filename
        print(f"PASS FASTA {filename}: records={len(names)}, aligned_width={width}")

    tree_text = (HERE / "gorilla_acres.treefile").read_text()
    tree_ids = re.findall(r"ACRE_[0-9]+", tree_text)
    assert len(tree_ids) == EXPECTED_TIPS
    assert len(set(tree_ids)) == EXPECTED_TIPS
    assert set(tree_ids) == expected
    iqtree_log = (HERE / "gorilla_acres.log").read_text()
    assert "BEST SCORE FOUND" in iqtree_log
    assert "Analysis results written to" in iqtree_log
    assert "Date and Time:" in iqtree_log
    print(f"PASS tree: unique_tips={len(set(tree_ids))}, IQ-TREE completion recorded")

    stems = [
        "tree_gorilla_clusters_radial",
        "tree_gorilla_clusters_rectangular",
        "tree_gorilla_chromosomes_radial",
        "tree_gorilla_chromosomes_rectangular",
    ]
    for stem in stems:
        for extension in ["png", "pdf"]:
            path = HERE / f"{stem}.{extension}"
            assert path.stat().st_size > 10_000, path
    print("PASS figures: four views x PNG/PDF")

    print("SHA256")
    for filename in [
        "gorilla_acres.treefile",
        "gorilla_tip_metadata.tsv",
        "gorilla_de_novo_clusters.tsv",
        "cluster_sweep_parameters.json",
    ]:
        print(f"{filename}\t{digest(HERE / filename)}")
    print("ALL VALIDATIONS PASSED")


if __name__ == "__main__":
    main()
