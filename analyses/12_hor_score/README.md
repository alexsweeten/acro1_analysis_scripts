# ACRE-scale HOR score

`scripts/compute_hor_score.py` is a generalized reconstruction of the scorer
used for the gorilla and HPRC results. Complete ACREs, rather than individual
ACRO1 monomers, are the repeating units in this adaptation.

```bash
python scripts/compute_hor_score.py \
  --metadata units.tsv \
  --fasta units.fasta \
  --output-dir results \
  --id-column acre_id \
  --array-column array_id \
  --order-column unit_in_array \
  --cluster-column cluster \
  --identity-thresholds 90,92,96 \
  --min-block-length 4 \
  --workers 8 \
  --seed 20261105
```

Outputs include per-unit scores, per-array block summaries, cluster summaries,
and a run manifest. `compute_gorilla_hor_score.py` is retained as the original
gorilla-specific implementation. The archived compressed HPRC score table is
included for verification because the original human generator was not
retained.

