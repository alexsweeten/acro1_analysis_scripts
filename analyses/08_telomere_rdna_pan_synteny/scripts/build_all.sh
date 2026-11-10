#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace_dir="$(cd "$script_dir/.." && pwd)"
renderer="$workspace_dir/synteny_chr13_primate_comparisons_10mb_compressed/render_pairwise_comparisons.R"

if (( $# > 0 )); then
  chromosomes=("$@")
else
  chromosomes=(hsa13 hsa14 hsa15 hsa21 hsa22)
fi

source_root_for() {
  case "$1" in
    hsa13) echo "$workspace_dir/synteny_chr13_primate_comparisons_10mb_compressed" ;;
    hsa14|hsa15|hsa21|hsa22) echo "$workspace_dir/synteny_hsa14_15_21_22_primate_comparisons_10mb_compressed/$1" ;;
    *) echo "Unknown chromosome: $1" >&2; return 1 ;;
  esac
}

for hsa in "${chromosomes[@]}"; do
  root="$script_dir/$hsa"
  source_root="$(source_root_for "$hsa")"
  mkdir -p "$root/sequences" "$root/masked_sequences" "$root/raw_alignments" \
    "$root/alignments" "$root/figures"

  Rscript "$script_dir/prepare_masked_inputs.R" "$source_root" "$root"

  while IFS=$'\t' read -r species_id label assembly sequence fasta bed endpoint cluster_start cluster_end; do
    if [[ "$species_id" == "species_id" ]]; then
      continue
    fi
    source_fasta="$script_dir/$fasta"
    samtools faidx -@ 2 -o "$root/sequences/$species_id.fa" \
      "$source_fasta" "$sequence:1-$endpoint"
    bedtools maskfasta -fi "$root/sequences/$species_id.fa" \
      -bed "$root/masks/$species_id.censat_mask.bed" \
      -fo "$root/masked_sequences/$species_id.masked.fa" -mc N
  done < "$root/species_config.tsv"

  for query_id in gorilla chimpanzee bonobo sumatran_orangutan bornean_orangutan; do
    minimap2 -t 4 -x asm5 -c --cs=short --secondary=yes -N 5 \
      -o "$root/raw_alignments/human_target_vs_$query_id.paf" \
      "$root/masked_sequences/human.masked.fa" \
      "$root/masked_sequences/$query_id.masked.fa"
    minimap2 -t 4 -x asm5 -c --cs=short --secondary=yes -N 5 \
      -o "$root/raw_alignments/${query_id}_target_vs_human.paf" \
      "$root/masked_sequences/$query_id.masked.fa" \
      "$root/masked_sequences/human.masked.fa"
  done

  Rscript "$script_dir/filter_reciprocal_alignments.R" "$root"
  Rscript "$renderer" "$root" "$hsa" acro_only unique_reciprocal
done
