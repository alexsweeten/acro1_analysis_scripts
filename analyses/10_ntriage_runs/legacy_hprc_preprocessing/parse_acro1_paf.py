#!/usr/bin/env python3
"""Filter sensitive ACRO1 mappings into monomer anchors for N-TRIAGE."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Hit:
    target: str
    strand: str
    target_start: int
    target_end: int
    query_start: int
    query_end: int
    matches: int
    alignment_length: int
    mapq: int
    alignment_score: int
    cigar: str

    @property
    def query_span(self) -> int:
        return self.query_end - self.query_start

    @property
    def target_span(self) -> int:
        return self.target_end - self.target_start

    @property
    def identity(self) -> float:
        return self.matches / self.alignment_length if self.alignment_length else 0.0


def parse_tags(fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for value in fields[12:]:
        pieces = value.split(":", 2)
        if len(pieces) == 3:
            tags[pieces[0]] = pieces[2]
    return tags


def overlaps_same_locus(left: Hit, right: Hit) -> bool:
    overlap = min(left.target_end, right.target_end) - max(left.target_start, right.target_start)
    return overlap > 0 and overlap / min(left.target_span, right.target_span) >= 0.50


def inferred_full_start(hit: Hit) -> int:
    query_offset = hit.query_start if hit.strand == "+" else 147 - hit.query_end
    return max(0, hit.target_start - query_offset)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--min-query-span", type=int, default=120)
    parser.add_argument("--min-target-span", type=int, default=120)
    parser.add_argument("--max-target-span", type=int, default=155)
    parser.add_argument("--min-identity", type=float, default=0.72)
    parser.add_argument("--min-score", type=int, default=80)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    targets: set[str] = set()
    with args.bed.open() as handle:
        for raw in handle:
            if raw.strip():
                targets.add(raw.split("\t", 1)[0])
    if len(targets) != 2794:
        raise ValueError(f"Expected 2,794 BED targets; observed {len(targets):,}")

    raw_hits: list[Hit] = []
    rejection: Counter[str] = Counter()
    total_alignments = 0
    with args.paf.open() as handle:
        for line_number, raw in enumerate(handle, start=1):
            if not raw.strip():
                continue
            total_alignments += 1
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"PAF line {line_number} has fewer than 12 fields")
            query, query_length = fields[0], int(fields[1])
            target = fields[5]
            if query != "ACRO1" or query_length != 147:
                rejection["wrong_query"] += 1
                continue
            if target not in targets:
                rejection["unknown_target"] += 1
                continue
            tags = parse_tags(fields)
            hit = Hit(
                target=target,
                strand=fields[4],
                target_start=int(fields[7]),
                target_end=int(fields[8]),
                query_start=int(fields[2]),
                query_end=int(fields[3]),
                matches=int(fields[9]),
                alignment_length=int(fields[10]),
                mapq=int(fields[11]),
                alignment_score=int(tags.get("AS", fields[9])),
                cigar=tags.get("cg", ""),
            )
            reason = ""
            if hit.query_span < args.min_query_span:
                reason = "short_query_span"
            elif not (args.min_target_span <= hit.target_span <= args.max_target_span):
                reason = "target_span_outside_range"
            elif hit.identity < args.min_identity:
                reason = "low_identity"
            elif hit.alignment_score < args.min_score:
                reason = "low_alignment_score"
            if reason:
                rejection[reason] += 1
            else:
                raw_hits.append(hit)

    ranked = sorted(
        raw_hits,
        key=lambda h: (-h.alignment_score, -h.identity, -h.query_span, h.target, h.strand, h.target_start),
    )
    kept: list[Hit] = []
    by_target_strand: dict[tuple[str, str], list[Hit]] = defaultdict(list)
    for hit in ranked:
        locus = by_target_strand[(hit.target, hit.strand)]
        if any(
            overlaps_same_locus(hit, prior)
            or abs(inferred_full_start(hit) - inferred_full_start(prior)) <= 6
            for prior in locus
        ):
            rejection["overlapping_alternative"] += 1
            continue
        locus.append(hit)
        kept.append(hit)
    kept.sort(key=lambda h: (h.target, inferred_full_start(h), h.strand))

    bed_path = args.out_dir / "acro1_monomer_anchors.bed"
    detail_path = args.out_dir / "acro1_monomer_anchor_details.tsv.gz"
    import gzip
    with bed_path.open("w") as bed, gzip.open(detail_path, "wt", newline="") as detail:
        columns = [
            "sample", "start", "end", "name", "score", "strand", "aligned_start",
            "aligned_end", "inferred_start", "query_start", "query_end", "query_span",
            "target_span", "matches", "alignment_length", "identity", "mapq",
            "alignment_score", "cigar",
        ]
        writer = csv.DictWriter(detail, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for number, hit in enumerate(kept, start=1):
            start = inferred_full_start(hit)
            end = start + 147
            name = f"{hit.target}_monomer_{number}"
            bed.write(f"{hit.target}\t{start}\t{end}\t{name}\t{hit.alignment_score}\t{hit.strand}\n")
            writer.writerow({
                "sample": hit.target, "start": start, "end": end, "name": name,
                "score": hit.alignment_score, "strand": hit.strand,
                "aligned_start": hit.target_start, "aligned_end": hit.target_end,
                "inferred_start": start, "query_start": hit.query_start,
                "query_end": hit.query_end, "query_span": hit.query_span,
                "target_span": hit.target_span, "matches": hit.matches,
                "alignment_length": hit.alignment_length, "identity": f"{hit.identity:.6f}",
                "mapq": hit.mapq, "alignment_score": hit.alignment_score, "cigar": hit.cigar,
            })

    exact_pairs = 0
    adjacent_pairs = 0
    for hits in by_target_strand.values():
        starts = sorted(inferred_full_start(hit) for hit in hits)
        adjacent_pairs += max(0, len(starts) - 1)
        exact_pairs += sum(b - a in (144, 147) for a, b in zip(starts, starts[1:]))
    qc_rows: list[tuple[str, object]] = [
        ("paf_alignments", total_alignments),
        ("passing_before_overlap_resolution", len(raw_hits)),
        ("retained_anchor_loci", len(kept)),
        ("arrays_with_anchor", len({hit.target for hit in kept})),
        ("retained_plus", sum(hit.strand == "+" for hit in kept)),
        ("retained_minus", sum(hit.strand == "-" for hit in kept)),
        ("adjacent_same_array_strand_pairs", adjacent_pairs),
        ("exact_144_or_147_pairs", exact_pairs),
        ("exact_pair_fraction", f"{exact_pairs / adjacent_pairs:.6f}" if adjacent_pairs else "NA"),
    ]
    qc_rows.extend((f"rejected_{reason}", count) for reason, count in sorted(rejection.items()))
    with (args.out_dir / "acro1_monomer_anchor_qc.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(qc_rows)

    print(f"Retained {len(kept):,} nonredundant ACRO1 anchors in {len({h.target for h in kept}):,} arrays")
    print(f"Exact 144/147-bp adjacent-pair fraction: {exact_pairs / adjacent_pairs:.6f}")


if __name__ == "__main__":
    main()
