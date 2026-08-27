"""scanpy: the per-group mean summaries an analyst gets from (i) sc.get.obs_df(...).groupby(...).mean(),
(ii) sc.get.aggregate(func='mean') and (iii) sc.pl.DotPlot(...).dot_color_df on normalize_total(1e4)+log1p data,
checked against the paper's lognorm_mean re-implementation (benchmark_aggregation.m_lognorm), and scored."""
import time, numpy as np, pandas as pd, hdf5plugin, anndata as ad, scipy.sparse as sp, scanpy as sc
from score_profiles import OUTDIR, BASE, relabel, score, PANEL
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_aggregation import m_lognorm, m_cpm, m_raw
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
RES = OUTDIR / "results"; RES.mkdir(exist_ok=True)
a = ad.read_h5ad(BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad")
a.obs["cell_type"] = pd.Categorical(relabel(a.obs["cell_type"].values))
lab = a.obs["cell_type"].astype(str).values; types = sorted(pd.unique(lab)); genes = a.var_names.values
X = sp.csr_matrix(a.X).astype(np.float32); expressed = np.asarray(X.sum(0)).ravel() > 0
P_ref = m_lognorm(X, lab, types)                                   # paper's re-implementation
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
log(f"scanpy {sc.__version__}; data log1p(CP10k) ready")

# (ii) sc.get.aggregate
agg = sc.get.aggregate(a, by="cell_type", func="mean")
P_agg = np.vstack([np.asarray(agg[t].layers["mean"]).ravel() for t in types])
log(f"sc.get.aggregate(func='mean') vs paper lognorm_mean: max|diff| = {np.abs(P_agg - P_ref).max():.3e}")

# (i) sc.get.obs_df + groupby, in gene blocks (a full dense frame would be 12k x 59k)
P_obs = np.zeros_like(P_ref); B = 4000
for j in range(0, a.n_vars, B):
    df = sc.get.obs_df(a, keys=list(genes[j:j+B]) + ["cell_type"])
    g = df.groupby("cell_type", observed=True).mean()
    P_obs[:, j:j+B] = g.loc[types, genes[j:j+B]].values
log(f"sc.get.obs_df().groupby().mean() vs paper lognorm_mean: max|diff| = {np.abs(P_obs - P_ref).max():.3e}")

# (iii) sc.pl.DotPlot on the 22 held-out markers (+ fraction expressing)
panel = [g for g in PANEL["HCC"] if g in set(genes)]
dp = sc.pl.DotPlot(a, var_names=panel, groupby="cell_type")
dc = dp.dot_color_df.loc[types, panel].values; ds = dp.dot_size_df.loc[types, panel].values
gidx = {g: i for i, g in enumerate(genes)}; ref22 = P_ref[:, [gidx[g] for g in panel]]
log(f"sc.pl.DotPlot.dot_color_df vs paper lognorm_mean (22 markers): max|diff| = {np.abs(dc - ref22).max():.3e}")
dot_arg = [types[i] for i in dc.argmax(0)]; size_arg = [types[i] for i in ds.argmax(0)]
tab = pd.DataFrame(dict(gene=panel, target=[PANEL["HCC"][g] for g in panel], dotplot_color_argmax=dot_arg,
                        dotplot_fraction_argmax=size_arg))
tab["color_correct"] = (tab.dotplot_color_argmax == tab.target).astype(int); tab["fraction_correct"] = (tab.dotplot_fraction_argmax == tab.target).astype(int)
tab.to_csv(RES / "scanpy_dotplot_markers.csv", index=False)
log(f"dotplot colour (mean log) argmax correct {tab.color_correct.sum()}/{len(tab)}; dot size (fraction expressing) argmax correct {tab.fraction_correct.sum()}/{len(tab)}")

rows = []
for name, what, P in [("sc.get.aggregate_mean", "sc.pp.normalize_total(1e4)+log1p; sc.get.aggregate(by='cell_type', func='mean')", P_agg),
                      ("obs_df_groupby_mean", "sc.pp.normalize_total(1e4)+log1p; sc.get.obs_df().groupby('cell_type').mean()", P_obs),
                      ("paper_lognorm_mean", "benchmark_aggregation.m_lognorm (re-implementation)", P_ref),
                      ("paper_cpm_mean", "benchmark_aggregation.m_cpm (re-implementation)", m_cpm(X, lab, types)),
                      ("paper_raw_mean", "benchmark_aggregation.m_raw (re-implementation)", m_raw(X, lab, types))]:
    r, mt = score(P, types, genes, expressed); rows.append(dict(tool="scanpy", path=name, computes=what, **r))
    log(f"{name:25} top {r['top_type']:>20} {r['top_share']:5.1%}  markers {r['markers_correct']}/{r['markers_total']} perm p={r['perm_p']:.4f}")
R = pd.DataFrame(rows); R["max_abs_diff_vs_paper_lognorm"] = [np.abs(P_agg-P_ref).max(), np.abs(P_obs-P_ref).max(), 0, np.nan, np.nan]
R.to_csv(RES / "scanpy_results.csv", index=False); log("wrote scanpy_results.csv")
