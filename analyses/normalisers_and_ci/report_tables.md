### Table 1. All methods, six real matrices (published values for the original six; this run for the new ten)

| method | top argmax share, median (range) | argmax in highest-coverage type | ρ(coverage, share) median; p<0.05 | held-out markers | perm p<0.05 | argmax identical to |
|---|---|---|---|---|---|---|
| raw mean | 73% (66%–98%) | 4/6 | +0.89; 5/6 | 92/198 | 6/6 | — |
| z-scaled mean (scanpy `scale`, clip 10) | 73% (65%–98%) | 4/6 | +0.86; 5/6 | 93/198 | 6/6 | — |
| within-cell rank mean | 66% (40%–99%) | 5/6 | +0.88; 6/6 | 93/198 | 5/6 | none |
| within-cell rank / nnz mean | 43% (33%–98%) | 5/6 | +0.89; 5/6 | 100/198 | 5/6 | none |
| detection-rate mean (binarised) | 52% (33%–100%) | 6/6 | +0.89; 5/6 | 87/198 | 5/6 | none |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | 37% (28%–95%) | 5/6 | +0.86; 4/6 | 107/198 | 6/6 | none |
| log-normalised mean (scanpy / Seurat / ArchR) | 36% (28%–83%) | 4/6 | +0.70; 2/6 | 105/198 | 6/6 | — |
| TF-IDF mean, Signac method 3 | 36% (28%–83%) | 4/6 | +0.70; 2/6 | 105/198 | 6/6 | lognorm_mean |
| log1p(CPM) mean, then per-gene max-normalisation | 36% (28%–83%) | 4/6 | +0.70; 2/6 | 105/198 | 6/6 | lognorm_mean |
| Pearson-residual mean (θ = 100) | 27% (22%–35%) | 2/6 | +0.60; 1/6 | 110/198 | 5/6 | none |
| coverage-matched thinning, raw mean | 28% (17%–31%) | 3/6 | +0.75; 3/6 | 107/198 | 5/6 | — |
| pooled CPM (sum per type → CPM) | 25% (24%–32%) | 0/6 | +0.57; 1/6 | 117/198 | 6/6 | — |
| pooled sum → log-normalised | 25% (24%–32%) | 0/6 | +0.57; 1/6 | 117/198 | 6/6 | pseudobulk_cpm |
| **CPM mean** | 25% (19%–30%) | 3/6 | +0.50; 1/6 | 106/198 | 5/6 | — |
| median-total scaling mean | 25% (19%–30%) | 3/6 | +0.50; 1/6 | 106/198 | 5/6 | cpm_mean |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | 25% (19%–30%) | 3/6 | +0.50; 1/6 | 106/198 | 5/6 | cpm_mean |

### Table 2. Per-matrix top argmax share (top type) for every method

| method | HCC prox | HCC enha | COAD prox | COAD enha | NSCLC prox | NSCLC enha |
|---|---|---|---|---|---|---|
| raw mean | 67% Endothelial * | 76% Endothelial * | 95% Mast* | 98% Mast* | 66% B cell | 71% B cell |
| z-scaled mean (scanpy `scale`, clip 10) | 65% Endothelial * | 75% Endothelial * | 95% Mast* | 98% Mast* | 65% B cell | 70% B cell |
| within-cell rank mean | 74% Endothelial * | 47% Endothelial * | 97% Mast* | 99% Mast* | 59% B cell | 40% T NK* |
| within-cell rank / nnz mean | 46% Endothelial * | 33% Hepatocyte | 92% Mast* | 98% Mast* | 34% T NK* | 41% T NK* |
| detection-rate mean (binarised) | 55% Endothelial * | 33% Endothelial * | 95% Mast* | 100% Mast* | 41% T NK* | 49% T NK* |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | 34% Endothelial * | 28% Hepatocyte | 91% Mast* | 95% Mast* | 37% T NK* | 36% T NK* |
| log-normalised mean (scanpy / Seurat / ArchR) | 34% Hepatocyte | 28% Hepatocyte | 60% Mast* | 83% Mast* | 37% T NK* | 36% T NK* |
| TF-IDF mean, Signac method 3 | 34% Hepatocyte | 28% Hepatocyte | 60% Mast* | 83% Mast* | 37% T NK* | 36% T NK* |
| log1p(CPM) mean, then per-gene max-normalisation | 34% Hepatocyte | 28% Hepatocyte | 60% Mast* | 83% Mast* | 37% T NK* | 36% T NK* |
| Pearson-residual mean (θ = 100) | 31% Hepatocyte | 25% Hepatocyte | 22% B cell | 35% B cell | 28% T NK* | 27% T NK* |
| coverage-matched thinning, raw mean | 31% Hepatocyte | 29% NK cytotoxic | 17% Macrophage | 23% Mast* | 28% T NK* | 28% T NK* |
| pooled CPM (sum per type → CPM) | 32% Hepatocyte | 26% Hepatocyte | 24% Macrophage | 32% B cell | 24% B cell | 25% B cell |
| pooled sum → log-normalised | 32% Hepatocyte | 26% Hepatocyte | 24% Macrophage | 32% B cell | 24% B cell | 25% B cell |
| **CPM mean** | 30% Hepatocyte | 25% Hepatocyte | 19% Macrophage | 21% Mast* | 25% T NK* | 25% T NK* |
| median-total scaling mean | 30% Hepatocyte | 25% Hepatocyte | 19% Macrophage | 21% Mast* | 25% T NK* | 25% T NK* |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | 30% Hepatocyte | 25% Hepatocyte | 19% Macrophage | 21% Mast* | 25% T NK* | 25% T NK* |

\* = the highest-coverage type (highest median true per-cell coverage).

### Table 3. Per-matrix held-out marker concordance (permutation p) for every method

| method | HCC prox | HCC enha | COAD prox | COAD enha | NSCLC prox | NSCLC enha | total |
|---|---|---|---|---|---|---|---|
| raw mean | 15/22 (p=0.000) | 14/22 (p=0.000) | 9/41 (p=0.010) | 9/41 (p=0.001) | 23/36 (p=0.000) | 22/36 (p=0.000) | 92/198 |
| z-scaled mean (scanpy `scale`, clip 10) | 16/22 (p=0.000) | 14/22 (p=0.000) | 9/41 (p=0.009) | 9/41 (p=0.001) | 23/36 (p=0.000) | 22/36 (p=0.000) | 93/198 |
| within-cell rank mean | 15/22 (p=0.000) | 18/22 (p=0.000) | 8/41 (p=0.018) | 5/41 (p=1.000) | 24/36 (p=0.000) | 23/36 (p=0.000) | 93/198 |
| within-cell rank / nnz mean | 18/22 (p=0.000) | 20/22 (p=0.000) | 9/41 (p=0.011) | 5/41 (p=0.874) | 25/36 (p=0.000) | 23/36 (p=0.000) | 100/198 |
| detection-rate mean (binarised) | 18/22 (p=0.000) | 14/22 (p=0.000) | 8/41 (p=0.040) | 5/41 (p=1.000) | 23/36 (p=0.000) | 19/36 (p=0.000) | 87/198 |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | 20/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.017) | 10/41 (p=0.000) | 23/36 (p=0.000) | 24/36 (p=0.000) | 107/198 |
| log-normalised mean (scanpy / Seurat / ArchR) | 19/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.025) | 9/41 (p=0.018) | 23/36 (p=0.000) | 24/36 (p=0.000) | 105/198 |
| TF-IDF mean, Signac method 3 | 19/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.026) | 9/41 (p=0.022) | 23/36 (p=0.000) | 24/36 (p=0.000) | 105/198 |
| log1p(CPM) mean, then per-gene max-normalisation | 19/22 (p=0.000) | 21/22 (p=0.000) | 9/41 (p=0.028) | 9/41 (p=0.021) | 23/36 (p=0.000) | 24/36 (p=0.000) | 105/198 |
| Pearson-residual mean (θ = 100) | 20/22 (p=0.000) | 21/22 (p=0.000) | 12/41 (p=0.006) | 10/41 (p=0.052) | 23/36 (p=0.000) | 24/36 (p=0.000) | 110/198 |
| coverage-matched thinning, raw mean | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.030) | 9/41 (p=0.063) | 24/36 (p=0.000) | 24/36 (p=0.000) | 107/198 |
| pooled CPM (sum per type → CPM) | 20/22 (p=0.000) | 21/22 (p=0.000) | 16/41 (p=0.000) | 14/41 (p=0.001) | 23/36 (p=0.000) | 23/36 (p=0.000) | 117/198 |
| pooled sum → log-normalised | 20/22 (p=0.000) | 21/22 (p=0.000) | 16/41 (p=0.000) | 14/41 (p=0.001) | 23/36 (p=0.000) | 23/36 (p=0.000) | 117/198 |
| **CPM mean** | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.032) | 9/41 (p=0.081) | 23/36 (p=0.000) | 24/36 (p=0.000) | 106/198 |
| median-total scaling mean | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.032) | 9/41 (p=0.087) | 23/36 (p=0.000) | 24/36 (p=0.000) | 106/198 |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | 19/22 (p=0.000) | 21/22 (p=0.000) | 10/41 (p=0.029) | 9/41 (p=0.080) | 23/36 (p=0.000) | 24/36 (p=0.000) | 106/198 |

### Table 4. Cell-level bootstrap (B = 200, cells resampled with replacement within each type) for the original six methods

Top argmax share: point estimate, percentile 95 % interval, basic (bias-corrected) 95 % interval, bootstrap bias; held-out marker count: point estimate and percentile 95 % interval; paired difference vs CPM mean (same resamples).

| matrix | method | top share | percentile 95 % CI | basic 95 % CI | bias | top type (stability) | markers | 95 % CI | Δshare vs CPM [95 % CI] | Δmarkers vs CPM [95 % CI] |
|---|---|---|---|---|---|---|---|---|---|---|
| HCC count | raw mean | 67.1% | 55.8%–68.7% | 65.4%–78.3% | -4.3% | Endothelial stromal (100%) | 15/22 | 14–17 | +37.0% [+27.4%, +40.5%] | -4 [-5, -2] |
| HCC count | z-scaled mean (scanpy `scale`, clip 10) | 65.3% | 54.2%–66.8% | 63.9%–76.5% | -4.0% | Endothelial stromal (100%) | 16/22 | 14–17 | +35.2% [+26.1%, +38.9%] | -3 [-5, -2] |
| HCC count | log-normalised mean (scanpy / Seurat / ArchR) | 33.6% | 30.4%–32.4% | 34.8%–36.9% | -2.2% | Hepatocyte (100%) | 19/22 | 19–20 | +3.5% [+2.8%, +4.1%] | +0 [+0, +1] |
| HCC count | coverage-matched thinning, raw mean | 30.3% | 27.5%–29.2% | 31.5%–33.2% | -2.0% | Hepatocyte (100%) | 19/22 | 18–20 | +0.2% [-0.4%, +0.9%] | +0 [-1, +1] |
| HCC count | pooled CPM (sum per type → CPM) | 31.6% | 28.1%–30.6% | 32.7%–35.2% | -2.2% | Hepatocyte (100%) | 20/22 | 18–20 | +1.6% [+0.3%, +2.3%] | +1 [-1, +1] |
| HCC count | **CPM mean** | 30.1% | 27.3%–28.8% | 31.4%–32.8% | -2.0% | Hepatocyte (100%) | 19/22 | 19–20 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| HCC weighted | raw mean | 75.6% | 63.0%–81.9% | 69.3%–88.2% | -1.6% | Endothelial stromal (100%) | 14/22 | 13–16 | +50.4% [+38.4%, +57.6%] | -7 [-8, -4] |
| HCC weighted | z-scaled mean (scanpy `scale`, clip 10) | 75.1% | 63.0%–81.0% | 69.1%–87.1% | -1.2% | Endothelial stromal (100%) | 14/22 | 13–16 | +49.8% [+39.1%, +56.6%] | -7 [-8, -5] |
| HCC weighted | log-normalised mean (scanpy / Seurat / ArchR) | 27.7% | 26.3%–27.9% | 27.6%–29.2% | -0.6% | Hepatocyte (100%) | 21/22 | 20–21 | +2.5% [+2.0%, +3.6%] | +0 [-1, +1] |
| HCC weighted | coverage-matched thinning, raw mean | 29.3% | 26.0%–29.3% | 29.3%–32.6% | -1.6% | NK cytotoxic T (100%) | 21/22 | 20–21 | +4.1% [+1.4%, +5.0%] | +0 [-1, +0] |
| HCC weighted | pooled CPM (sum per type → CPM) | 26.4% | 24.8%–26.4% | 26.4%–28.0% | -0.8% | Hepatocyte (100%) | 21/22 | 19–21 | +1.2% [+0.5%, +1.8%] | +0 [-2, +1] |
| HCC weighted | **CPM mean** | 25.2% | 23.8%–25.0% | 25.4%–26.6% | -0.8% | Hepatocyte (100%) | 21/22 | 20–21 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| COAD count | raw mean | 94.8% | 91.0%–93.2% | 96.3%–98.6% | -2.5% | Mast (100%) | 9/41 | 9–10 | +76.0% [+73.5%, +75.8%] | -1 [-3, +4] |
| COAD count | z-scaled mean (scanpy `scale`, clip 10) | 94.9% | 91.3%–93.5% | 96.3%–98.5% | -2.3% | Mast (100%) | 9/41 | 9–10 | +76.1% [+74.0%, +76.1%] | -1 [-4, +4] |
| COAD count | log-normalised mean (scanpy / Seurat / ArchR) | 60.1% | 55.7%–57.5% | 62.7%–64.4% | -3.5% | Mast (100%) | 9/41 | 8–11 | +41.3% [+38.3%, +40.2%] | -1 [-4, +4] |
| COAD count | coverage-matched thinning, raw mean | 17.4% | 15.0%–15.8% | 19.0%–19.8% | -2.0% | Macrophage (100%) | 10/41 | 6–13 | -1.4% [-2.5%, -1.7%] | +0 [-4, +5] |
| COAD count | pooled CPM (sum per type → CPM) | 23.9% | 21.0%–22.3% | 25.5%–26.8% | -2.3% | Macrophage (90%) | 16/41 | 12–17 | +5.1% [+3.5%, +4.9%] | +6 [+2, +10] |
| COAD count | **CPM mean** | 18.8% | 17.1%–17.7% | 19.9%–20.5% | -1.4% | Macrophage (98%) | 10/41 | 6–13 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| COAD weighted | raw mean | 97.6% | 97.2%–97.8% | 97.4%–98.1% | -0.1% | Mast (100%) | 9/41 | 7–10 | +76.7% [+78.5%, +80.2%] | +0 [-3, +4] |
| COAD weighted | z-scaled mean (scanpy `scale`, clip 10) | 97.7% | 97.3%–97.9% | 97.5%–98.1% | -0.1% | Mast (100%) | 9/41 | 7–10 | +76.8% [+78.6%, +80.2%] | +0 [-3, +3] |
| COAD weighted | log-normalised mean (scanpy / Seurat / ArchR) | 82.8% | 80.7%–82.2% | 83.5%–84.9% | -1.3% | Mast (100%) | 9/41 | 8–10 | +61.9% [+62.2%, +64.3%] | +0 [-2, +4] |
| COAD weighted | coverage-matched thinning, raw mean | 22.8% | 19.2%–20.4% | 25.2%–26.4% | -3.0% | Mast (100%) | 11/41 | 6–13 | +1.9% [+1.0%, +2.2%] | +2 [-2, +5] |
| COAD weighted | pooled CPM (sum per type → CPM) | 31.7% | 28.7%–30.0% | 33.3%–34.6% | -2.3% | B cell (100%) | 14/41 | 12–18 | +10.7% [+10.2%, +12.0%] | +5 [+2, +10] |
| COAD weighted | **CPM mean** | 20.9% | 17.5%–18.9% | 23.0%–24.4% | -2.7% | Mast (99%) | 9/41 | 5–12 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| NSCLC count | raw mean | 66.1% | 59.9%–70.4% | 61.8%–72.3% | -0.8% | B cell (100%) | 23/36 | 22–24 | +40.8% [+36.0%, +46.3%] | +0 [-3, +1] |
| NSCLC count | z-scaled mean (scanpy `scale`, clip 10) | 64.6% | 58.4%–69.0% | 60.2%–70.9% | -0.8% | B cell (100%) | 23/36 | 21–24 | +39.4% [+34.3%, +45.0%] | +0 [-3, +1] |
| NSCLC count | log-normalised mean (scanpy / Seurat / ArchR) | 36.9% | 33.7%–38.1% | 35.7%–40.1% | -1.0% | T NK (100%) | 23/36 | 23–24 | +11.7% [+9.7%, +14.0%] | +0 [-2, +1] |
| NSCLC count | coverage-matched thinning, raw mean | 28.2% | 25.8%–27.6% | 28.8%–30.6% | -1.5% | T NK (100%) | 24/36 | 22–25 | +2.9% [+1.8%, +3.6%] | +1 [-2, +1] |
| NSCLC count | pooled CPM (sum per type → CPM) | 24.1% | 22.9%–23.9% | 24.4%–25.4% | -0.8% | B cell (98%) | 23/36 | 23–25 | -1.1% [-1.3%, +0.0%] | +0 [-2, +1] |
| NSCLC count | **CPM mean** | 25.3% | 23.5%–24.5% | 26.0%–27.0% | -1.3% | T NK (100%) | 23/36 | 23–25 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |
| NSCLC weighted | raw mean | 71.3% | 65.8%–74.4% | 68.1%–76.7% | -0.5% | B cell (100%) | 22/36 | 21–23 | +46.1% [+42.2%, +51.0%] | -2 [-4, +1] |
| NSCLC weighted | z-scaled mean (scanpy `scale`, clip 10) | 70.4% | 64.7%–73.9% | 66.8%–76.1% | -0.6% | B cell (100%) | 22/36 | 21–23 | +45.2% [+40.9%, +50.3%] | -2 [-4, +1] |
| NSCLC weighted | log-normalised mean (scanpy / Seurat / ArchR) | 35.5% | 32.5%–36.3% | 34.7%–38.5% | -1.0% | T NK (100%) | 24/36 | 23–25 | +10.4% [+9.3%, +12.9%] | +0 [-2, +2] |
| NSCLC weighted | coverage-matched thinning, raw mean | 27.7% | 25.3%–27.9% | 27.4%–30.0% | -1.2% | T NK (100%) | 25/36 | 23–26 | +2.6% [+1.8%, +4.5%] | +1 [-2, +3] |
| NSCLC weighted | pooled CPM (sum per type → CPM) | 24.6% | 23.3%–24.5% | 24.7%–25.8% | -0.7% | B cell (100%) | 23/36 | 22–25 | -0.5% [-0.3%, +1.3%] | -1 [-3, +1] |
| NSCLC weighted | **CPM mean** | 25.1% | 22.9%–24.0% | 26.3%–27.3% | -1.7% | T NK (100%) | 24/36 | 22–26 | +0.0% [+0.0%, +0.0%] | +0 [+0, +0] |

Thinning-only noise (10 independent binomial thinnings, no cell resampling), depth-matched method: top-share SD HCC count 0.0009 (markers 19-20), HCC weighted 0.0016 (markers 21-21), COAD count 0.0009 (markers 9-15), COAD weighted 0.0020 (markers 9-11), NSCLC count 0.0014 (markers 23-24), NSCLC weighted 0.0015 (markers 23-25).

### Table 5. Paired tests across the six matrices, each method vs CPM mean (point estimates; exact Wilcoxon signed-rank and exact sign test; n = 6 pairs, smallest attainable two-sided Wilcoxon p = 0.031, one-sided 0.016)

| method | metric | matrices with method > CPM | median Δ | Wilcoxon p (two-sided) | Wilcoxon p (method > CPM) | sign-test p |
|---|---|---|---|---|---|---|
| raw mean | top share | 6/6 | +0.483 | 0.031 | 0.016 | 0.031 |
| raw mean | marker fraction | 0/6 | -0.040 | 0.125 | 1.000 | 0.125 |
| z-scaled mean (scanpy `scale`, clip 10) | top share | 6/6 | +0.475 | 0.031 | 0.016 | 0.031 |
| z-scaled mean (scanpy `scale`, clip 10) | marker fraction | 0/6 | -0.040 | 0.125 | 1.000 | 0.125 |
| within-cell rank mean | top share | 6/6 | +0.386 | 0.031 | 0.016 | 0.031 |
| within-cell rank mean | marker fraction | 1/6 | -0.073 | 0.094 | 0.984 | 0.219 |
| within-cell rank / nnz mean | top share | 6/6 | +0.157 | 0.031 | 0.016 | 0.031 |
| within-cell rank / nnz mean | marker fraction | 1/6 | -0.037 | 0.312 | 0.891 | 0.219 |
| detection-rate mean (binarised) | top share | 6/6 | +0.247 | 0.031 | 0.016 | 0.031 |
| detection-rate mean (binarised) | marker fraction | 0/6 | -0.073 | 0.062 | 1.000 | 0.062 |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | top share | 6/6 | +0.115 | 0.031 | 0.016 | 0.031 |
| TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4)) | marker fraction | 2/6 | +0.000 | 0.750 | 0.375 | 1.000 |
| log-normalised mean (scanpy / Seurat / ArchR) | top share | 6/6 | +0.110 | 0.031 | 0.016 | 0.031 |
| log-normalised mean (scanpy / Seurat / ArchR) | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| TF-IDF mean, Signac method 3 | top share | 6/6 | +0.110 | 0.031 | 0.016 | 0.031 |
| TF-IDF mean, Signac method 3 | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| log1p(CPM) mean, then per-gene max-normalisation | top share | 6/6 | +0.110 | 0.031 | 0.016 | 0.031 |
| log1p(CPM) mean, then per-gene max-normalisation | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| Pearson-residual mean (θ = 100) | top share | 6/6 | +0.023 | 0.031 | 0.016 | 0.031 |
| Pearson-residual mean (θ = 100) | marker fraction | 3/6 | +0.012 | 0.250 | 0.125 | 0.250 |
| coverage-matched thinning, raw mean | top share | 5/6 | +0.023 | 0.094 | 0.047 | 0.219 |
| coverage-matched thinning, raw mean | marker fraction | 1/6 | +0.000 | 1.000 | 0.500 | 1.000 |
| pooled CPM (sum per type → CPM) | top share | 4/6 | +0.014 | 0.156 | 0.078 | 0.688 |
| pooled CPM (sum per type → CPM) | marker fraction | 3/6 | +0.023 | 0.250 | 0.125 | 0.625 |
| pooled sum → log-normalised | top share | 4/6 | +0.014 | 0.156 | 0.078 | 0.688 |
| pooled sum → log-normalised | marker fraction | 3/6 | +0.023 | 0.250 | 0.125 | 0.625 |
| median-total scaling mean | top share | 1/6 | +0.000 | 1.000 | 0.500 | 1.000 |
| median-total scaling mean | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | top share | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
| TF-IDF mean, Signac method 2 (TF·log(1+IDF)) | marker fraction | 0/6 | +0.000 | 1.000 | 1.000 | 1.000 |
