# scatac-aggregation-benchmark

Code for **"Linear-Scale Per-Cell Normalisation Mitigates Library-Size-Associated
Variation in Aggregated Single-Cell ATAC-seq Gene-Activity Profiles."**

Gene-activity matrices from single-cell ATAC-seq are used as lineage priors for
cluster annotation, marker weighting and reference construction. Each use begins by
aggregating cells into one profile per cell type, and these matrices hold per-cell
values whose totals differ between cell types by factors of two to forty. The scale
on which that aggregation is performed is rarely reported and has not been
evaluated. This repository compares six aggregation strategies on six gene-activity
matrices from three tumour datasets (liver, colorectal, lung), on a deliberately
mislabelled negative control, on a matched RNA–ATAC sample, and on a controlled
series in which the difference in per-cell totals was imposed from one- to
forty-fold.

## The six strategies, on the HCC count matrix

| strategy | per-type profile | top share | markers |
|---|---|---:|---:|
| raw mean | mean of unnormalised values | 67.1% | 15/22 |
| z-scaled mean | `scanpy.pp.scale` per gene, then mean | 65.3% | 16/22 |
| log-normalised mean | scale each cell to 10⁴, log(1+x), mean | 33.6% | 19/22 |
| pooled CPM | sum per type, scale to 10⁶ | 31.6% | 20/22 |
| **CPM mean** | scale each cell to 10⁴, mean on the linear scale | **30.1%** | 19/22 |
| coverage-matched | binomial thinning to a common row sum, then raw mean | 30.5% | 19/22 |

"Top share" is the share of genes whose per-gene maximum falls in the single most
frequent cell type — with T types and no library-size effect it should sit near 1/T.

Across all six matrices, averaging unnormalised or per-gene z-scaled values places a
median 73% of genome-wide per-gene maxima in a single cell type, against 25% for the
CPM mean; leave-one-sample-out ranges for the two do not overlap in any dataset.
Averaging log-transformed normalised values gives 36% and remains sensitive to the
imposed difference, reaching 83% in the dataset with the widest spread, whereas
linear-scale profiles stay between 22% and 30% across that range. Of 198 held-out
markers the raw mean places 92 in their own lineage, the z-scaled mean 93, the
log-normalised mean 105, the CPM mean 106 and pooled CPM 117.

The mechanism: a cell contributes 0 to the mean of a gene it does not detect, and
the per-type detection rate rises with coverage. Only a transform linear in the
per-cell fraction x/d cancels that term. The mean of log(1+CPM) instead tracks each
type's detection rate (ρ = 0.95 with coverage in the colorectal matrix).

## Setup

```bash
pip install -r requirements.txt        # or see environment.txt for exact versions
```

No data is committed. Download the deposited files listed in
[`data/README.md`](data/README.md) — checksums in
[`data/checksums.md`](data/checksums.md) — then point the scripts at them:

```bash
export SCRNA_ATAC_BASE=/path/to/your/data/root   # default: the repository root
export GENCODE_GTF=/path/to/gencode.v38.annotation.gtf
export GSE227265_FRAGMENTS=/path/to/GSE227265_fragments_AllSamples.tsv.gz
```

`data/README.md` gives the expected layout under `$SCRNA_ATAC_BASE` and says which
of the three audits needs which files — none needs the whole tree.

## Running

Build the gene-activity matrices from the deposited files:

```bash
python3 scripts/90_enhancer_aware_gene_activity.py   # HCC
python3 scripts/92_coad_gene_activity_rebuild.py     # colorectal
python3 scripts/94_luad_prior_and_reference.py       # NSCLC
```

Then the benchmark and its diagnostics:

```bash
python3 tools/benchmark_aggregation.py        # -> results/benchmark_aggregation.csv, _markers.csv, _depth_by_type.csv
python3 tools/check_hcc_labels.py             # -> results/hcc_label_check.csv
python3 tools/benchmark_wilcoxon_markers.py   # -> results/benchmark_wilcoxon.csv
python3 tools/revision_analyses.py            # -> results/revision/
python3 tools/make_benchmark_figure.py        # -> figures/ (gitignored)
```

Marker panels (`PANEL`), the matrix list (`MATRICES`), the permutation null and the
random seed (`SEED = 42`) are all defined in `tools/benchmark_aggregation.py`; every
other script imports them, so the strategy definitions are shared.

## Supporting analyses

Five independent analyses in [`analyses/`](analyses/INDEX.md), each with its own
scripts, summary tables and a `REPORT.md`: an imposed coverage-spread series, a
confirmation against the real Seurat / Signac / scanpy / SnapATAC2 code paths, a
matched RNA–ATAC ground truth, ten further normalisers with bootstrap intervals,
and a survey of what published pipelines actually do.

## Layout

```
tools/       the benchmark, the label check, the audits and the figure scripts
scripts/     gene-activity matrix construction from the deposited files
results/     the summary tables the benchmark writes
analyses/    five supporting analyses, each with scripts, tables and a REPORT.md
data/        accessions, expected layout and checksums (no data files are stored here)
```

## Naming

The manuscript's final terminology postdates most of this code. Identifiers in
scripts, CSV columns and directory names keep the older forms:

| in the code | in the manuscript |
|---|---|
| `depth` | coverage |
| `pseudobulk_cpm` | pooled CPM |
| `depth_matched` | coverage-matched |
| `proximal` | count matrix |
| `enhancer` | weighted matrix |
| `LUAD` | NSCLC |

## Superseded files

`tools/audit_*.py` and `results/marker_concordance*.csv`, `coad_depth_audit.csv` and
`zscale_*.csv` belong to an earlier version of this work and are kept for the
record. `marker_concordance.csv` was computed with the two HCC cluster labels that
`check_hcc_labels.py` shows to be wrong; files ending `_RAW_superseded.csv` are the
pre-correction versions.

Large derived matrices — the gene × cell-type profile tables written by
`analyses/real_tools/` — are not committed; re-run the scripts in that directory to
regenerate them.

## Licence

MIT — see [`LICENSE`](LICENSE). Author: Yen-Jung Chiu.
