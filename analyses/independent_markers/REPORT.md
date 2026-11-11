# Marker concordance with markers from independent references

**Question.** The held-out marker panels of the benchmark were chosen by the author. Does the
ordering of the aggregation strategies hold with markers taken from independent references?

**Markers.** CellTypist logistic-regression models trained on independent atlases
(celltypist 1.7.1; models as distributed by CellTypist):

| dataset | model | source atlas |
|---|---|---|
| PBMC | Immune_All_Low | Domínguez Conde et al., Science 2022 |
| brain | Adult_Human_PrefrontalCortex | Ma et al., Science 2022 |
| HCC | Healthy_Human_Liver | Guilliams et al., Cell 2022 |
| COAD | Human_Colorectal_Cancer | Lee et al., Nat Genet 2020 |
| NSCLC | Human_Lung_Atlas | Sikkema et al., Nat Med 2023 |

Each of our cell types is mapped to model classes by fixed name rules written before any result
was seen (`marker_mapping.csv`). Specificity score of gene g for our type t: mean coefficient over
t's classes minus the largest mean coefficient over any other type. Panel: the top 10 (and 20)
genes with a positive score. Excluded: every gene used to assign the labels (PBMC and brain
labelling panels, the 72-gene COAD and 35-gene NSCLC scoring sets, every gene of the HCC label
check) and genes absent from either gene-activity matrix. Scoring and permutation null (10,000
shuffles) as for the held-out panels. Panels: `marker_panels.csv`.

**Overlap with the author's panels** (top 10): PBMC 0/80, brain 0/50, HCC 4/60, COAD 13/80,
NSCLC 2/60.

**Result** (`independent_markers.csv`, top 10, ten matrices; null ≈ 10 per matrix):

| strategy | all five datasets | without COAD |
|---|---|---|
| raw mean | 207 / 660 | 177 / 500 |
| CP10K mean | 252 / 660 | 230 / 500 |
| log-normalized mean | 261 / 660 | 231 / 500 |
| pooled CPM | 263 / 660 | 231 / 500 |

Largest differences: brain (raw 15 vs 26–30 of 50 in either matrix) and HCC (30–31 vs 37–40 of
60). COAD is the exception: CP10K mean 10 and 12 of 80 (at the null) vs 15 for the raw mean;
COAD labels were named on the unnormalized scale. With 20 markers per type the raw mean is lowest
or tied everywhere except COAD and, by one marker, the log-normalized mean in the PBMC weighted
matrix.

Run: `python3 01_independent_markers.py` (needs the CellTypist models in `~/.celltypist`).
