# Non-tumor tissues and multi-donor RNA references

The six aggregation strategies of `tools/benchmark_aggregation.py`, imported unchanged, applied to

| dataset | what | units |
|---|---|---|
| GSE184462 (adult samples) | human chromatin accessibility atlas, 72 fragment files | 26 non-tumor tissues, 487,097 cells |
| 10x `human_brain_3k` | brain multiome, labels from RNA only | 1 donor, 2,901 cells, 5 types |
| GSE194122 (multiome) | bone-marrow multiome, gene-activity matrix as deposited by its authors | 10 donors, 69,249 cells, 22 types |

Run order (paths at the top of each script point to the raw downloads):

```bash
python3 01_build_catlas.py --list                 # tissues with fragment files on disk
python3 01_build_catlas.py liver lung ...         # count + distance-weighted gene activity per tissue
python3 02_benchmark_catlas.py                    # -> catlas_benchmark.csv, catlas_loso.csv, catlas_markers.csv, catlas_depth_by_type.csv
python3 03_build_brain3k.py                       # gene activity for the brain multiome
python3 04_label_brain_rna.py                     # RNA-only labels -> brain3k_rna_labels.csv
python3 05_benchmark_brain.py                     # -> brain3k_*.csv
python3 08_brain_share_by_type.py                 # -> brain3k_share_by_type.csv (share of maxima per type)
python3 06_benchmark_bmmc.py                      # -> bmmc_benchmark.csv, bmmc_lodo.csv, bmmc_depth_by_type.csv
python3 07_figure_normal_tissues.py               # figure + fig_normal_tissues_data.csv
```

Notes
- The two `UMB4540_snATAC_frontal_cortex` files of GSE184462 are the atlas samples `Human_brain_1` / `Human_brain_2`.
- Cell types with fewer than 100 cells in a tissue are dropped; the atlas labels are the depositors' own.
- The bone-marrow gene activity is deposited as log1p(count * s / total) with a batch-specific constant s; integer
  counts are recovered exactly (every entry within 1e-4 of an integer, asserted in `06_benchmark_bmmc.py`).
