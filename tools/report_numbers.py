"""Recompute and print every quantity the analysis reports, from the result
tables in results/ and analyses/.

This replaces the manuscript-checking script used locally: it derives the values
from the tables alone, so the repository does not need the manuscript source.
Run after any re-analysis and compare the printed values with the ones you cite.
"""
import math
from pathlib import Path
import pandas as pd
from scipy.stats import spearmanr

R = Path(__file__).resolve().parent.parent
res, rev, adds = R / "results", R / "results" / "revision", R / "analyses"
pct = lambda v: math.floor(v * 100 + 0.5)

B = pd.read_csv(res / "benchmark_aggregation.csv"); real = B[B.model != "misannotated"]
D = pd.read_csv(res / "benchmark_depth_by_type.csv")
M = pd.read_csv(res / "benchmark_markers.csv")
C = pd.read_csv(res / "benchmark_mechanism_coad.csv").sort_values("det")
NAME = {"LUAD": "NSCLC"}

print("== datasets (count matrices) ==")
for coh in ("HCC", "LUAD", "COAD"):
    g = D[(D.cohort == coh) & (D.model == "proximal")].sort_values("median_true_depth")
    lo, hi = g.iloc[0], g.iloc[-1]
    sh = {m: pct(B[(B.cohort == coh) & (B.model == "proximal") & (B.method == m)].top_share.iloc[0])
          for m in ("raw_mean", "lognorm_mean", "cpm_mean")}
    print(f"  {NAME.get(coh, coh):6} {int(g.n_cells.sum()):>7,} cells | {lo.cell_type} {lo.median_true_depth:.0f}"
          f" -> {hi.cell_type} {hi.median_true_depth:.0f} ({hi.median_true_depth/lo.median_true_depth:.1f}x)"
          f" | detection {g.detection_rate.min()*100:.1f}-{g.detection_rate.max()*100:.1f}%"
          f" | top share raw {sh['raw_mean']} log {sh['lognorm_mean']} CPM {sh['cpm_mean']}")

print("\n== top share and marker concordance, six matrices ==")
for m in ("raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "cpm_mean"):
    s = real[real.method == m]
    print(f"  {m:15} median {pct(s.top_share.median())}% ({pct(s.top_share.min())}-{pct(s.top_share.max())})"
          f" | modal type is the highest-coverage one in {int((s.top_type == s.deepest_type).sum())}/6"
          f" | rho median {s.rho_depth_share.median():+.2f} | markers {int(s.markers_correct.sum())}/198")
dm = {c: pct(B[(B.cohort == c) & (B.model == "proximal") & (B.method == "depth_matched")].top_share.iloc[0])
      for c in ("HCC", "COAD", "LUAD")}
print(f"  coverage-matched control (count matrices): HCC {dm['HCC']}%, COAD {dm['COAD']}%, NSCLC {dm['LUAD']}%")
print(f"  smallest permutation p in the table: {B.perm_p.min():.5g} (zeros: {int((B.perm_p == 0).sum())})")

print("\n== leave-one-sample-out top share ==")
L = pd.read_csv(rev / "loso_top_share.csv"); g = L.groupby(["cohort", "method"]).top_share.agg(["min", "max", "count"])
for c in ("COAD", "HCC", "LUAD"):
    for m in ("raw_mean", "lognorm_mean", "cpm_mean"):
        r = g.loc[(c, m)]
        print(f"  {NAME.get(c, c):6} {m:14} {r['min']*100:.1f}-{r['max']*100:.1f}% (n={int(r['count'])})")
    print(f"    raw and CPM ranges disjoint: {g.loc[(c,'raw_mean'),'min'] > g.loc[(c,'cpm_mean'),'max']}")

print("\n== mechanism, colorectal count matrix ==")
cd = D[(D.cohort == "COAD") & (D.model == "proximal")]
Cm = C.merge(cd[["cell_type", "median_true_depth"]], left_on="type", right_on="cell_type")
print(f"  detection rate {cd.detection_rate.min()*100:.1f}-{cd.detection_rate.max()*100:.1f}%,"
      f" rho with coverage {spearmanr(cd.median_true_depth, cd.detection_rate)[0]:+.2f}")
print(f"  per-type mean of log(1+CPM): rho with coverage {spearmanr(Cm.median_true_depth, Cm.log)[0]:+.2f},"
      f" {C.log.max()/C.log.min():.1f}-fold range")
print(f"  per-type mean of CPM: constant at {C.cpm.iloc[0]:.3f} across {len(C)} types")

print("\n== matrix integrality ==")
I = pd.read_csv(rev / "matrix_integrality_and_shares.csv")
for r in I.sort_values(["cohort", "model"]).itertuples():
    print(f"  {NAME.get(r.cohort, r.cohort):6} {r.model:9} integer={r.integer!s:5} "
          f"non-zero entries below 0.5: {r.frac_below_0_5:.3f}")

print("\n== HCC label sensitivity ==")
H = pd.read_csv(rev / "hcc_label_sensitivity.csv")
p = H.pivot_table(index=["model", "method"], columns="labels", values="markers_correct")
print(p.to_string())

print("\n== controlled depth series ==")
S = pd.read_csv(adds / "depth_spread" / "depth_spread_summary.csv")
for c in ("HCC", "LUAD"):
    for m in ("raw_mean", "lognorm_mean", "cpm_mean", "pseudobulk_cpm"):
        r = S[(S.cohort == c) & (S.method == m)].sort_values("spread")
        print(f"  {NAME.get(c, c):6} {m:15} " + "  ".join(
            f"{int(x.spread)}x:{pct(x.top_share_mean)}%({int(x.deep_wins)}/6)" for x in r.itertuples()))

print("\n== matched RNA-ATAC ==")
U = pd.read_csv(adds / "multiome" / "multiome_benchmark.csv")
T = pd.read_csv(adds / "multiome" / "multiome_rna_truth_share.csv")
print(f"  RNA reference, share in the modal type: "
      f"{T[(T.model=='proximal')&(T.cell_type=='pDC')].rna_argmax_share_all.iloc[0]*100:.1f}%")
for r in U.itertuples():
    print(f"  {r.model:9} {r.method:15} top share {r.top_share*100:5.1f}% | RNA agreement "
          f"{r.argmax_agree_informative:.3f} | per-gene rho {r.gene_rho_mean_informative:.3f} | "
          f"markers {r.markers_correct}/{r.markers_total} | annotation {r.annot_acc_pearson:.4f}")
