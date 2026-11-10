#!/usr/bin/env python3
"""Plot the HSA19-orthologous ACRO1 locus across primates.

Synteny ribbons are drawn only from the ten independent unique-sequence
anchors. Repeat-derived matches are displayed inside the inferred orthologous
locus but are never used to connect species.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Polygon, Rectangle
from matplotlib.lines import Line2D
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
INTERMEDIATE = ROOT / "intermediate"
FIGURES = ROOT / "figures"

ALL_ORDER = [
    "Human",
    "Chimpanzee",
    "Bonobo",
    "Gorilla",
    "Sumatran_orangutan",
    "Bornean_orangutan",
    "Siamang_hap1",
    "Siamang_hap2",
    "Hoolock_gibbon",
    "Nomascus_gibbon",
    "Rhesus_macaque",
    "Marmoset",
    "Golden_lion_tamarin",
    "Squirrel_monkey",
]

REPRESENTATIVE_ORDER = [
    "Human",
    "Chimpanzee",
    "Gorilla",
    "Sumatran_orangutan",
    "Siamang_hap1",
    "Nomascus_gibbon",
    "Rhesus_macaque",
    "Marmoset",
    "Squirrel_monkey",
]

DISPLAY_NAMES = {
    "Human": "Human",
    "Chimpanzee": "Chimpanzee",
    "Bonobo": "Bonobo",
    "Gorilla": "Gorilla",
    "Sumatran_orangutan": "Sumatran orangutan",
    "Bornean_orangutan": "Bornean orangutan",
    "Siamang_hap1": "Siamang, haplotype 1",
    "Siamang_hap2": "Siamang, haplotype 2",
    "Hoolock_gibbon": "Hoolock gibbon",
    "Nomascus_gibbon": "Nomascus gibbon",
    "Rhesus_macaque": "Rhesus macaque",
    "Marmoset": "Common marmoset",
    "Golden_lion_tamarin": "Golden lion tamarin",
    "Squirrel_monkey": "Squirrel monkey",
}

CLADE_COLORS = {
    "great_ape": "#274C77",
    "gibbon": "#6A4C93",
    "old_world_monkey": "#2A9D8F",
    "new_world_monkey": "#C65D21",
}

ANCHOR_COLORS = {
    "L500": "#0B5D4B",
    "L250": "#1F8068",
    "L100": "#45A98B",
    "L40": "#84CDB7",
    "L20": "#BFE8DC",
    "R20": "#FAD8B5",
    "R40": "#F5B477",
    "R100": "#ED8A45",
    "R250": "#D85A2A",
    "R500": "#A93226",
}

ACRO_COLOR = "#7349A4"
COMP_COLOR = "#E69F00"
MOBILE_COLOR = "#7C8798"
LOCUS_FILL = "#F7F7F7"


def read_inputs():
    locus = pd.read_csv(RESULTS / "orthologous_locus_matrix.tsv", sep="\t")
    architecture = pd.read_csv(RESULTS / "orthologous_locus_architecture.tsv", sep="\t")
    windows = pd.read_csv(RESULTS / "ortholog_window_manifest.tsv", sep="\t")
    anchors = pd.read_csv(RESULTS / "all_anchor_hits.tsv", sep="\t")
    hmm = pd.read_csv(RESULTS / "ACRO1_HMM_hits_at_orthologs.tsv", sep="\t")
    locus = locus.set_index("assembly_key", drop=False)
    architecture = architecture.set_index("assembly_key", drop=False)
    windows = windows.set_index("assembly_key", drop=False)
    anchors = anchors[anchors["selected_chain"].astype(str).str.lower() == "true"].copy()
    return locus, architecture, windows, anchors, hmm


def oriented_interval(start: float, end: float, locus_row) -> tuple[float, float]:
    center = (float(locus_row.locus_start) + float(locus_row.locus_end)) / 2.0
    if locus_row.orientation_to_human == "+":
        return start - center, end - center
    return center - end, center - start


def merge_intervals(intervals, max_gap=0):
    cleaned = sorted((float(a), float(b)) for a, b in intervals if b > a)
    if not cleaned:
        return []
    merged = [list(cleaned[0])]
    for start, end in cleaned[1:]:
        if start <= merged[-1][1] + max_gap:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(interval) for interval in merged]


def parse_component_paf(path: Path):
    hits = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip().split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            tstart = int(fields[7])
            tend = int(fields[8])
            nmatch = int(fields[9])
            alnlen = int(fields[10])
            coverage = (qend - qstart) / qlen if qlen else 0
            identity = nmatch / alnlen if alnlen else 0
            if coverage < 0.45 or identity < 0.60 or alnlen < 100:
                continue
            hits.append({"qname": qname, "tstart": tstart, "tend": tend})
    return hits


def collect_plot_segments(locus, architecture, windows, anchors, hmm):
    plot_anchors = defaultdict(dict)
    excluded_anchors = []
    segment_rows = []

    for _, hit in anchors.iterrows():
        key = hit.assembly_key
        if key not in locus.index:
            continue
        x0, x1 = oriented_interval(hit.assembly_start, hit.assembly_end, locus.loc[key])
        # Nomascus L20 is a known locally transposed hit 12.8 Mb away. It was
        # not used to bracket the inferred ortholog and is excluded from the
        # common plotting window.
        if max(abs(x0), abs(x1)) > 1_000_000:
            excluded_anchors.append(
                {
                    "assembly_key": key,
                    "anchor": hit.anchor,
                    "reason": "outside_common_plot_range_likely_local_transposition",
                    "oriented_start_bp": x0,
                    "oriented_end_bp": x1,
                }
            )
            continue
        plot_anchors[key][hit.anchor] = (x0, x1)
        segment_rows.append(
            {
                "assembly_key": key,
                "feature_type": "unique_anchor",
                "feature_name": hit.anchor,
                "oriented_start_bp": x0,
                "oriented_end_bp": x1,
            }
        )

    acro_segments = defaultdict(list)
    for key, group in hmm.groupby("assembly_key"):
        if key not in locus.index or key not in windows.index:
            continue
        win = windows.loc[key]
        central_start = float(locus.loc[key].locus_start)
        central_end = float(locus.loc[key].locus_end)
        raw_intervals = []
        for _, hit in group.iterrows():
            start = float(win.window_start) + min(float(hit.start), float(hit.end))
            end = float(win.window_start) + max(float(hit.start), float(hit.end))
            start = max(start, central_start)
            end = min(end, central_end)
            if end > start:
                raw_intervals.append(oriented_interval(start, end, locus.loc[key]))
        acro_segments[key] = merge_intervals(raw_intervals, max_gap=175)
        for x0, x1 in acro_segments[key]:
            segment_rows.append(
                {
                    "assembly_key": key,
                    "feature_type": "ACRO1_HMM_tract",
                    "feature_name": "ACRO1",
                    "oriented_start_bp": x0,
                    "oriented_end_bp": x1,
                }
            )

    component_segments = defaultdict(lambda: defaultdict(list))
    for key in locus.index:
        path = INTERMEDIATE / "per_window_component_mappings" / f"{key}.components.paf"
        if not path.exists() or key not in windows.index:
            continue
        win = windows.loc[key]
        central_start = float(locus.loc[key].locus_start)
        central_end = float(locus.loc[key].locus_end)
        grouped = defaultdict(list)
        for hit in parse_component_paf(path):
            if hit["qname"].startswith("ACRO1_"):
                continue
            category = "mobile_element" if hit["qname"].startswith(("L1", "MER")) else "ACRE_component"
            start = max(float(win.window_start) + hit["tstart"], central_start)
            end = min(float(win.window_start) + hit["tend"], central_end)
            if end > start:
                grouped[category].append(oriented_interval(start, end, locus.loc[key]))
        for category, intervals in grouped.items():
            component_segments[key][category] = merge_intervals(intervals, max_gap=250)
            for x0, x1 in component_segments[key][category]:
                segment_rows.append(
                    {
                        "assembly_key": key,
                        "feature_type": category,
                        "feature_name": category,
                        "oriented_start_bp": x0,
                        "oriented_end_bp": x1,
                    }
                )

    pd.DataFrame(segment_rows).to_csv(RESULTS / "HSA19_synteny_plot_segments.tsv", sep="\t", index=False)
    pd.DataFrame(excluded_anchors).to_csv(
        RESULTS / "HSA19_synteny_plot_excluded_anchors.tsv", sep="\t", index=False
    )
    return plot_anchors, acro_segments, component_segments, excluded_anchors


def add_interval(ax, x0, x1, y0, height, color, edgecolor="none", linewidth=0.0, zorder=3, alpha=1.0):
    min_width = 900.0
    center = (x0 + x1) / 2.0
    width = max(x1 - x0, min_width)
    ax.add_patch(
        Rectangle(
            (center - width / 2.0, y0),
            width,
            height,
            facecolor=color,
            edgecolor=edgecolor,
            linewidth=linewidth,
            zorder=zorder,
            alpha=alpha,
        )
    )


def draw_figure(order, basename, title_suffix, locus, architecture, plot_anchors, acro_segments, component_segments):
    n = len(order)
    height = 0.72 * n + 3.4
    fig, ax = plt.subplots(figsize=(17.2, height))
    y_by_key = {key: n - 1 - index for index, key in enumerate(order)}

    # Ribbons join the same unique anchor in adjacent plotted assemblies only.
    for upper, lower in zip(order[:-1], order[1:]):
        y_upper = y_by_key[upper]
        y_lower = y_by_key[lower]
        for anchor in ANCHOR_COLORS:
            if anchor not in plot_anchors[upper] or anchor not in plot_anchors[lower]:
                continue
            a0, a1 = plot_anchors[upper][anchor]
            b0, b1 = plot_anchors[lower][anchor]
            polygon = Polygon(
                [(a0, y_upper - 0.10), (a1, y_upper - 0.10),
                 (b1, y_lower + 0.10), (b0, y_lower + 0.10)],
                closed=True,
                facecolor=ANCHOR_COLORS[anchor],
                edgecolor="none",
                alpha=0.11,
                zorder=0,
            )
            ax.add_patch(polygon)

    for key in order:
        y = y_by_key[key]
        loc = locus.loc[key]
        arch = architecture.loc[key]
        span = float(loc.locus_span_bp)
        half = span / 2.0

        anchor_values = list(plot_anchors[key].values())
        line_start = min((x[0] for x in anchor_values), default=-half)
        line_end = max((x[1] for x in anchor_values), default=half)
        ax.plot([line_start, line_end], [y, y], color="#AEB7C2", linewidth=1.1, zorder=1)

        # Inferred orthologous interval, with a strong black outline.
        ax.add_patch(
            Rectangle(
                (-half, y - 0.25), span, 0.50,
                facecolor=LOCUS_FILL, edgecolor="#111111", linewidth=1.45, zorder=2,
            )
        )

        for anchor, (x0, x1) in plot_anchors[key].items():
            add_interval(ax, x0, x1, y - 0.15, 0.30, ANCHOR_COLORS[anchor], zorder=4)

        for x0, x1 in component_segments[key].get("mobile_element", []):
            add_interval(ax, x0, x1, y - 0.20, 0.105, MOBILE_COLOR, zorder=5)
        for x0, x1 in component_segments[key].get("ACRE_component", []):
            add_interval(ax, x0, x1, y - 0.075, 0.12, COMP_COLOR, zorder=6)
        for x0, x1 in acro_segments[key]:
            add_interval(ax, x0, x1, y + 0.07, 0.13, ACRO_COLOR, zorder=7)

        clade_color = CLADE_COLORS.get(str(loc.clade), "#333333")
        ax.scatter([-714_000], [y], s=33, color=clade_color, clip_on=False, zorder=8)

        # Right-side aligned text table.
        state = str(arch.orthologous_locus_state).replace("_", " ")
        ax.text(
            1.012, y, f"{span/1000:6.1f}", transform=ax.get_yaxis_transform(),
            ha="right", va="center", fontsize=8.3, color="#333333", clip_on=False,
        )
        ax.text(
            1.070, y, f"{float(arch.ACRO1_bp_in_inferred_locus)/1000:5.1f}",
            transform=ax.get_yaxis_transform(), ha="right", va="center",
            fontsize=8.3, color=ACRO_COLOR, clip_on=False,
        )
        ax.text(
            1.085, y, state, transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=8.3, color="#333333", clip_on=False,
        )

    # Anchor labels on the human row.
    human_y = y_by_key[order[0]]
    for anchor in ANCHOR_COLORS:
        if anchor not in plot_anchors[order[0]]:
            continue
        x0, x1 = plot_anchors[order[0]][anchor]
        ax.text(
            (x0 + x1) / 2, human_y + 0.33, anchor,
            ha="center", va="bottom", fontsize=7.2, color=ANCHOR_COLORS[anchor],
            rotation=45 if abs((x0 + x1) / 2) < 130_000 else 0,
        )

    ax.set_xlim(-720_000, 720_000)
    ax.set_ylim(-1.00, n + 0.50)
    ax.set_yticks([y_by_key[key] for key in order])
    ax.set_yticklabels([DISPLAY_NAMES[key] for key in order], fontsize=9.2)
    for label, key in zip(ax.get_yticklabels(), order):
        label.set_color(CLADE_COLORS.get(str(locus.loc[key].clade), "#222222"))
        if key == "Human":
            label.set_fontweight("bold")

    ticks = list(range(-600_000, 600_001, 200_000))
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{tick/1000:.0f}" for tick in ticks])
    ax.set_xlabel("Distance from the inferred locus center (kb; all tracks oriented to human)", fontsize=10)
    ax.axvline(0, color="#303030", linewidth=0.6, linestyle=(0, (2, 3)), zorder=0)
    ax.grid(axis="x", color="#E5E7EB", linewidth=0.6, zorder=-2)
    ax.tick_params(axis="y", length=0, pad=10)
    ax.tick_params(axis="x", labelsize=8.5)

    # Clade separators.
    for upper, lower in zip(order[:-1], order[1:]):
        if locus.loc[upper].clade != locus.loc[lower].clade:
            y_sep = (y_by_key[upper] + y_by_key[lower]) / 2.0
            ax.axhline(y_sep, color="#D1D5DB", linewidth=0.8, zorder=-1)

    ax.text(1.012, 1.008, "Locus\n(kb)", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.2, weight="bold")
    ax.text(1.070, 1.008, "ACRO1\n(kb)", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.2, weight="bold", color=ACRO_COLOR)
    ax.text(1.085, 1.008, "Locus state", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.2, weight="bold")

    title = "HSA19-orthologous ACRO1 locus across primates"
    if title_suffix:
        title += f" — {title_suffix}"
    ax.set_title(title, loc="left", fontsize=15.5, weight="bold", pad=35)
    ax.text(
        0, 1.025,
        "Unique-sequence anchors define synteny; repeat-derived features are shown only within the inferred locus",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=9.4, color="#4B5563",
    )

    legend = [
        Line2D([0], [0], color="#9CA3AF", lw=1.5, label="Unique-flank backbone"),
        Patch(facecolor="#4D9C88", alpha=0.30, edgecolor="none", label="Unique-anchor synteny ribbon"),
        Patch(facecolor=LOCUS_FILL, edgecolor="#111111", linewidth=1.2, label="Inferred orthologous interval"),
        Patch(facecolor=ACRO_COLOR, label="ACRO1 profile-HMM tract"),
        Patch(facecolor=COMP_COLOR, label="Non-ACRO ACRE component"),
        Patch(facecolor=MOBILE_COLOR, label="Associated mobile-element component"),
    ]
    ax.legend(
        handles=legend, loc="upper left", bbox_to_anchor=(0.0, -0.105), ncol=3,
        frameon=False, fontsize=8.4, handlelength=1.8, columnspacing=1.5,
    )
    ax.text(
        0.0, -0.205,
        "Common linear scale; locus lengths are exact. Nomascus L20 transposed hit (+12.8 Mb) omitted from plotting.",
        transform=ax.transAxes, ha="left", va="top", fontsize=7.7, color="#6B7280",
    )

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#9CA3AF")
    fig.subplots_adjust(left=0.17, right=0.76, top=0.88, bottom=0.22)

    for extension in ["png", "pdf", "svg"]:
        output = FIGURES / f"{basename}.{extension}"
        fig.savefig(output, dpi=320 if extension == "png" else None, facecolor="white")
    plt.close(fig)


def main():
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.titlecolor": "#16181D",
            "text.color": "#16181D",
        }
    )
    FIGURES.mkdir(parents=True, exist_ok=True)
    locus, architecture, windows, anchors, hmm = read_inputs()
    plot_anchors, acro_segments, component_segments, excluded = collect_plot_segments(
        locus, architecture, windows, anchors, hmm
    )
    draw_figure(
        ALL_ORDER,
        "HSA19_ACRO1_synteny_all14",
        "all 14 assemblies",
        locus,
        architecture,
        plot_anchors,
        acro_segments,
        component_segments,
    )
    draw_figure(
        REPRESENTATIVE_ORDER,
        "HSA19_ACRO1_synteny_representative9",
        "representative panel",
        locus,
        architecture,
        plot_anchors,
        acro_segments,
        component_segments,
    )
    print(f"Wrote HSA19 synteny figures to {FIGURES}")
    print(f"Plotted {sum(len(v) for v in plot_anchors.values())} unique-anchor blocks")
    print(f"Excluded {len(excluded)} distal/transposed anchor hit(s)")


if __name__ == "__main__":
    main()
