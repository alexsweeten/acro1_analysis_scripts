# Analysis directories

Each numbered directory corresponds to a Methods analysis. Scripts are copied
from the final analysis branch rather than from exploratory figure-layout
iterations. Compact TSV files are included when they record fixed inputs,
manual decisions, or the values reported in the manuscript.

The code falls into three categories:

1. **Parameterized analysis:** accepts command-line arguments or documented
   environment variables and can be rerun after acquiring the public inputs.
2. **Workflow provenance:** records the exact analysis logic, but requires
   reconstructing the original multi-stage input directory from public data.
3. **Manual audit:** includes both code and the frozen row-level decisions.

| Directory | Category | Important note |
| --- | --- | --- |
| `01_acre_extraction_and_alleles` | Workflow provenance | AniAnn's and RepeatMasker are external dependencies. |
| `02_species_balanced_phylogeny` | Parameterized analysis | Set `ACREAGE_ROOT` and `HSA19_ORTHOLOGY_ROOT`. |
| `03_gorilla_phylogeny_and_clusters` | Workflow provenance | Includes the all-1,232-unit tree and G1--G6 scripts. |
| `04_orf_analysis` | Workflow provenance | Large tract FASTAs are not redistributed. |
| `05_polr2a_composition` | Parameterized analysis | Compact allele and POLR2A inputs are bundled. |
| `06_prdm9_motif` | Parameterized analysis | Set input paths with CLI options or `ACRO1_DATA_ROOT`. |
| `07_acro_hsat1a` | Parameterized analysis | Run `Rscript run_six_species_analysis.R --help`. |
| `08_telomere_rdna_pan_synteny` | Workflow provenance | Requires assembly-matched CenSat BED and genome FASTA files. |
| `09_hsa19_orthology` | Workflow provenance | Large ortholog windows and PAFs are omitted. |
| `10_ntriage_runs` | Mixed | Requires NTRiage v0.3.0 and the HPRC inputs. |
| `11_hprc_transitions` | Manual audit | Publish the `audit/` tables with the code. |
| `12_hor_score` | Parameterized analysis | Generic reconstructed scorer plus archived result summaries. |
| `13_divergence` | Mixed | HPRC Mash workflow is directly parameterized. |
| `14_promoter_conservation` | Mixed | PRO-seq source provenance remains unresolved. |
| `15_assembly_tables` | Parameterized utility | Assembly and HPRC accession tables are bundled. |

