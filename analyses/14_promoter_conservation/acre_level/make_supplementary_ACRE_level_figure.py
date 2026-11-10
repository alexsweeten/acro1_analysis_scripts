#!/usr/bin/env python3
"""Build the supplementary cluster-specific ACRE promoter figure."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parent
SOURCE_SCRIPT = ROOT / "figure6_plotting_base.py"
OUT = ROOT / "figures"


def load_source_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("figure6_acre_source", SOURCE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {SOURCE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    source = load_source_module()
    source.configure_style()
    OUT.mkdir(parents=True, exist_ok=True)

    profile = pd.read_csv(source.PROFILE, sep="\t")
    values = pd.read_csv(source.GROUP_VALUES, sep="\t")
    summary = pd.read_csv(source.SUMMARY, sep="\t")
    one_sample = pd.read_csv(source.ONE_SAMPLE, sep="\t")
    pairwise = pd.read_csv(source.PAIRWISE, sep="\t")

    figure = plt.figure(figsize=(10.8, 18.6))
    grid = figure.add_gridspec(
        3,
        1,
        height_ratios=[1.04, 1.04, 0.82],
        hspace=0.40,
    )
    axis_profile = figure.add_subplot(grid[0])
    axis_effect = figure.add_subplot(grid[1])
    axis_scaled = figure.add_subplot(grid[2])

    source.localized_profile_panel(axis_profile, profile)
    source.violin_panel(axis_effect, values, summary, one_sample)
    axis_effect.set_ylim(-3.55, 3.85)
    source.normalized_panel(axis_scaled, summary, pairwise)

    for label, axis in zip(("A", "B", "C"), (axis_profile, axis_effect, axis_scaled)):
        axis.text(
            -0.105,
            1.08,
            label,
            transform=axis.transAxes,
            fontsize=28,
            fontweight="bold",
            va="top",
        )

    figure.suptitle(
        "Cluster-specific L1PA10 promoter conservation among H1-H3 ACREs",
        fontsize=24,
        fontweight="bold",
        y=0.988,
    )
    figure.subplots_adjust(left=0.12, right=0.975, top=0.947, bottom=0.055)

    pdf = OUT / "SupplementaryFigure_ACRE_level_L1PA10_promoter_conservation.pdf"
    png = OUT / "SupplementaryFigure_ACRE_level_L1PA10_promoter_conservation.png"
    figure.savefig(pdf, bbox_inches="tight")
    figure.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    plt.close(figure)
    print(pdf)
    print(png)


if __name__ == "__main__":
    main()
