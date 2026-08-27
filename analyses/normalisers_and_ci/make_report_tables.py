"""Assemble the markdown tables for REPORT.md from the merged CSVs (prints to stdout and
writes report_tables.md). The interpretation in REPORT.md is written by hand."""
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent
A = pd.read_csv(HERE / "benchmark_all_methods.csv")
C = pd.read_csv(HERE / "benchmark_bootstrap_ci.csv")
T = pd.read_csv(HERE / "benchmark_paired_tests.csv")
PUB = pd.read_csv(HERE.parents[2] / "results" / "benchmark_aggregation.csv"); PUB = PUB[PUB.model != "misannotated"]
ORDER = ["raw_mean", "zscale_mean", "rank_mean", "rankfrac_mean", "detection_mean", "tfidf1_mean", "lognorm_mean",
         "tfidf3_mean", "lognorm_maxnorm", "pearson_mean", "depth_matched", "pseudobulk_cpm", "pseudobulk_logcpm",
         "cpm_mean", "medscale_mean", "tfidf2_mean"]
ORIG = {"raw_mean", "zscale_mean", "pseudobulk_cpm", "depth_matched", "cpm_mean", "lognorm_mean"}
LABEL = {"raw_mean": "raw mean", "zscale_mean": "z-scale mean (scanpy `scale`, clip 10)", "rank_mean": "within-cell rank mean",
         "rankfrac_mean": "within-cell rank / nnz mean", "detection_mean": "detection-rate mean (binarised)",
         "tfidf1_mean": "TF-IDF mean, Signac method 1 (log1p(TF·IDF·1e4))",
         "lognorm_mean": "log1p(CPM) mean (scanpy / Seurat / ArchR)", "tfidf3_mean": "TF-IDF mean, Signac method 3",
         "lognorm_maxnorm": "log1p(CPM) mean, then per-gene max-normalisation", "pearson_mean": "Pearson-residual mean (θ = 100)",
         "depth_matched": "depth-matched thinning, raw mean", "pseudobulk_cpm": "pseudobulk sum → CPM",
         "pseudobulk_logcpm": "pseudobulk sum → log1p(CPM)", "cpm_mean": "**CPM mean (ours)**",
         "medscale_mean": "median-total scaling mean", "tfidf2_mean": "TF-IDF mean, Signac method 2 (TF·log(1+IDF))"}
MATS = [("HCC", "proximal"), ("HCC", "enhancer"), ("COAD", "proximal"), ("COAD", "enhancer"), ("LUAD", "proximal"), ("LUAD", "enhancer")]
out = []

def src_rows(m):
    """Original six: the published rows (results/benchmark_aggregation.csv); new: this run."""
    return PUB[PUB.method == m] if m in ORIG else A[A.method == m]

out.append("### Table 1. All methods, six real matrices (published values for the original six; this run for the new ten)\n")
out.append("| method | top argmax share, median (range) | argmax in deepest type | ρ(depth, share) median; p<0.05 | held-out markers | perm p<0.05 | argmax identical to |")
out.append("|---|---|---|---|---|---|---|")
for m in ORDER:
    s = src_rows(m)
    ident = sorted(set(A[A.method == m].argmax_identical_to.dropna())) if m not in ORIG else []
    ident = [i for i in ident if i != "none"]
    out.append(f"| {LABEL[m]} | {s.top_share.median():.0%} ({s.top_share.min():.0%}–{s.top_share.max():.0%}) | "
               f"{(s.top_type == s.deepest_type).sum()}/{len(s)} | {s.rho_depth_share.median():+.2f}; {(s.p_depth_share < 0.05).sum()}/{len(s)} | "
               f"{s.markers_correct.sum()}/{s.markers_total.sum()} | {(s.perm_p < 0.05).sum()}/{len(s)} | "
               f"{', '.join(ident) if ident else ('—' if m in ORIG else 'none')} |")

out.append("\n### Table 2. Per-matrix top argmax share (top type) for every method\n")
out.append("| method | " + " | ".join(f"{c} {mo[:4]}" for c, mo in MATS) + " |")
out.append("|---|" + "---|" * len(MATS))
for m in ORDER:
    s = src_rows(m).set_index(["cohort", "model"])
    cells = []
    for c, mo in MATS:
        r = s.loc[(c, mo)]
        cells.append(f"{r.top_share:.0%} {r.top_type.replace('_', ' ')[:12]}{'*' if r.top_type == r.deepest_type else ''}")
    out.append(f"| {LABEL[m]} | " + " | ".join(cells) + " |")
out.append("\n\\* = the deepest type (highest median true per-cell depth).")

out.append("\n### Table 3. Per-matrix held-out marker concordance (permutation p) for every method\n")
out.append("| method | " + " | ".join(f"{c} {mo[:4]}" for c, mo in MATS) + " | total |")
out.append("|---|" + "---|" * (len(MATS) + 1))
for m in ORDER:
    s = src_rows(m).set_index(["cohort", "model"])
    cells = [f"{int(s.loc[(c, mo)].markers_correct)}/{int(s.loc[(c, mo)].markers_total)} (p={s.loc[(c, mo)].perm_p:.3f})" for c, mo in MATS]
    out.append(f"| {LABEL[m]} | " + " | ".join(cells) + f" | {int(s.markers_correct.sum())}/{int(s.markers_total.sum())} |")

out.append("\n### Table 4. Cell-level bootstrap (B = 200, cells resampled with replacement within each type) for the original six methods\n")
out.append("Top argmax share: point estimate, percentile 95 % interval, basic (bias-corrected) 95 % interval, bootstrap bias; "
           "held-out marker count: point estimate and percentile 95 % interval; paired difference vs CPM mean (same resamples).\n")
out.append("| matrix | method | top share | percentile 95 % CI | basic 95 % CI | bias | top type (stability) | markers | 95 % CI | Δshare vs CPM [95 % CI] | Δmarkers vs CPM [95 % CI] |")
out.append("|---|---|---|---|---|---|---|---|---|---|---|")
for c, mo in MATS:
    for m in ["raw_mean", "zscale_mean", "lognorm_mean", "depth_matched", "pseudobulk_cpm", "cpm_mean"]:
        r = C[(C.cohort == c) & (C.model == mo) & (C.method == m)]
        if r.empty: continue
        r = r.iloc[0]
        out.append(f"| {c} {mo} | {LABEL[m]} | {r.top_share:.1%} | {r.top_share_ci_lo:.1%}–{r.top_share_ci_hi:.1%} | "
                   f"{r.top_share_basic_lo:.1%}–{r.top_share_basic_hi:.1%} | {r.top_share_boot_bias:+.1%} | "
                   f"{r.top_type.replace('_', ' ')} ({r.top_type_stability:.0%}) | {int(r.markers_correct)}/{int(r.markers_total)} | "
                   f"{r.markers_ci_lo:.0f}–{r.markers_ci_hi:.0f} | {r.diff_share_vs_cpm:+.1%} [{r.diff_share_ci_lo:+.1%}, {r.diff_share_ci_hi:+.1%}] | "
                   f"{int(r.diff_markers_vs_cpm):+d} [{r.diff_markers_ci_lo:+.0f}, {r.diff_markers_ci_hi:+.0f}] |")
dm = C[C.method == "depth_matched"]
out.append("\nThinning-only noise (10 independent binomial thinnings, no cell resampling), depth-matched method: top-share SD "
           + ", ".join(f"{r.cohort} {r.model} {r.thinning_sd:.4f} (markers {r.thinning_markers_range})" for r in dm.itertuples()) + ".")

out.append("\n### Table 5. Paired tests across the six matrices, each method vs CPM mean (point estimates; exact Wilcoxon signed-rank and exact sign test; n = 6 pairs, smallest attainable two-sided Wilcoxon p = 0.031, one-sided 0.016)\n")
out.append("| method | metric | matrices with method > CPM | median Δ | Wilcoxon p (two-sided) | Wilcoxon p (method > CPM) | sign-test p |")
out.append("|---|---|---|---|---|---|---|")
for m in ORDER:
    for metric in ["top_share", "marker_fraction"]:
        r = T[(T.method == m) & (T.metric == metric)]
        if r.empty: continue
        r = r.iloc[0]
        out.append(f"| {LABEL[m]} | {metric.replace('_', ' ')} | {int(r.n_positive)}/{int(r.n_pairs)} | {r.median_diff:+.3f} | "
                   f"{r.wilcoxon_p_two_sided:.3f} | {r.wilcoxon_p_greater:.3f} | {r.sign_test_p:.3f} |")
txt = "\n".join(out)
(HERE / "report_tables.md").write_text(txt + "\n")
print(txt)
