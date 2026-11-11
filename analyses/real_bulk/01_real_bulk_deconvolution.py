"""
01: Deconvolution of real bulk ATAC-seq with known composition (added 2026-11-11 in response to
the editor's comment on practical significance).

Bulk: five PBMC samples of the EPIC-ATAC study (Gabriel et al., eLife 2024; Zenodo
10.5281/zenodo.12571947, CC BY 4.0): PBMC_counts.txt (peak x sample counts, hg38) and
PBMC_cell_fractions.txt (flow-cytometry proportions). Downloaded unchanged to
data/raw/epic_atac_pbmc_bulk/.
Bulk gene activity: the counts of all peaks overlapping a gene's window (gene body + 2 kb
upstream, GENCODE v38), i.e. the window of the count model, for the genes of the PBMC count matrix.
References: the 10x PBMC multiome count matrix with RNA-derived labels (as in the paper), one
reference per aggregation strategy (functions imported from tools/benchmark_aggregation.py).
Signature genes: 200 per type by specificity ratio on the CP10K-mean reference, as in the
deconvolution simulation, restricted to genes with bulk counts; shared by every strategy.
Fit: NNLS of each bulk sample (CPM over signature genes) on each reference; coefficients rescaled
to sum to one; estimates merged to the flow-cytometry groups: B; monocytes (CD14 + CD16);
CD4 T (naive + non-naive + Treg in flow); CD8 T; NK; DC (DC + pDC).
Out: real_bulk_estimates.csv, real_bulk_metrics.csv, run.log
"""
import os, sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
from scipy.optimize import nnls
warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from benchmark_aggregation import m_raw, m_zscale, m_pseudobulk_cpm, m_depth_matched, m_cpm, m_lognorm  # noqa
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "/mnt/10t/scrna_atac"))
BULK = BASE / "data/raw/epic_atac_pbmc_bulk"
GTF = os.environ.get("GENCODE_GTF", "/mnt/10t/holiday/hnsc_analysis/gencode.v38.annotation.gtf")
REF = BASE / "data/processed/multiome_pbmc10k/pbmc10k_geneactivity_proximal.h5ad"
LAB = ROOT / "analyses" / "multiome" / "rna_labels.csv"
SEED, NSIG = 42, 200
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def P_of(r): return r[0] if isinstance(r, tuple) else r

GROUP = {"B": "B", "Mono_CD14": "Monocytes", "Mono_CD16": "Monocytes", "T_CD4": "CD4 T", "T_CD8": "CD8 T",
         "NK": "NK", "DC": "DC", "pDC": "DC"}
FLOW = {"Bcells": "B", "Monocytes": "Monocytes", "Naive_CD4_Tcells": "CD4 T", "Non_Naive_CD4_Tcells": "CD4 T",
        "Tregs": "CD4 T", "Naive_CD8_Tcells": "CD8 T", "Non_Naive_CD8_Tcells": "CD8 T", "NKcells": "NK", "Dendritic": "DC"}

# ── flow-cytometry truth ──
F = pd.read_csv(BULK / "PBMC_cell_fractions.txt", sep="\t", index_col=0)
F = F.drop(columns=[c for c in F.columns if c == "cell_type"])
assert set(F.index) == set(FLOW), sorted(F.index)
truth = F.groupby(F.index.map(FLOW)).sum()                 # groups x samples
log(f"flow truth (column sums {truth.sum().round(3).to_dict()}):\n{truth.round(3)}")

# ── bulk peak counts -> gene windows ──
B = pd.read_csv(BULK / "PBMC_counts.txt", sep="\t", index_col=0)
assert set(B.columns) == set(truth.columns), (B.columns, truth.columns)
B = B[truth.columns]                                        # align samples by name
loc = B.index.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
assert loc.notna().all().all()
pk = pd.DataFrame({"chrom": loc[0].values, "s": loc[1].astype(int).values, "e": loc[2].astype(int).values})
cnt = B.values.astype(np.float64)
a = ad.read_h5ad(REF)
lab_df = pd.read_csv(LAB, index_col=0); lab_df = lab_df[lab_df.cell_type != "dropped"]
a = a[lab_df.index.intersection(a.obs_names)].copy()
lab = lab_df.loc[a.obs_names, "cell_type"].values.astype(str)
types = sorted(set(lab)); assert set(types) == set(GROUP), types
genes = a.var_names.astype(str).values
win = {}
with open(GTF) as fh:
    for line in fh:
        if line[0] == "#": continue
        p = line.split("\t")
        if p[2] != "gene": continue
        name = p[8].split('gene_name "')[1].split('"')[0]
        if name in win: continue
        s, e = int(p[3]), int(p[4])
        win[name] = (p[0], s - 2000, e) if p[6] == "+" else (p[0], s, e + 2000)
Gb = np.zeros((len(genes), cnt.shape[1]))
by_chr = {c: g.sort_values("s") for c, g in pk.assign(i=np.arange(len(pk))).groupby("chrom")}
cums = {c: (g.s.values, g.e.values, np.vstack([np.zeros(cnt.shape[1]), np.cumsum(cnt[g.i.values], 0)])) for c, g in by_chr.items()}
assert all(np.all(np.diff(v[1]) >= 0) for v in cums.values()), "peak ends not sorted with starts"
n_win = 0
for gi, g in enumerate(genes):
    if g not in win: continue
    c, ws, we = win[g]
    if c not in cums: continue
    s, e, cs = cums[c]
    lo = np.searchsorted(e, ws, side="right"); hi = np.searchsorted(s, we, side="left")
    if hi > lo: Gb[gi] = cs[hi] - cs[lo]
    n_win += 1
log(f"bulk: {len(pk)} peaks, {cnt.shape[1]} samples; gene windows for {n_win}/{len(genes)} genes; "
    f"fraction of peak counts inside some window (sample 1): {Gb[:, 0].sum() / cnt[:, 0].sum():.2f} (windows overlap)")

# ── references ──
X = sp.csr_matrix(a.X).astype(np.float32)
rng = np.random.default_rng(SEED)
REFS = {"raw_mean": P_of(m_raw(X, lab, types)), "zscale_mean": P_of(m_zscale(X, lab, types)),
        "lognorm_mean": P_of(m_lognorm(X, lab, types)), "pseudobulk_cpm": P_of(m_pseudobulk_cpm(X, lab, types)),
        "cpm_mean": P_of(m_cpm(X, lab, types)), "depth_matched": P_of(m_depth_matched(X, lab, types, rng))}
C = REFS["cpm_mean"]; c0 = 0.01 * C[C > 0].mean()
has_bulk = Gb.sum(1) > 0
sig = set()
for k in range(len(types)):
    ratio = (C[k] + c0) / (np.delete(C, k, 0).max(0) + c0)
    ratio[~has_bulk] = -np.inf
    sig.update(np.argsort(-ratio)[:NSIG].tolist())
sig = np.array(sorted(sig))
log(f"signature genes: {len(sig)}")
Y = Gb[sig] / Gb[sig].sum(0) * 1e6                          # bulk CPM over signature genes

est_rows, met_rows = [], []
groups = list(truth.index)
for name, P in REFS.items():
    A = P[:, sig].T.astype(np.float64)                     # genes x types
    E = []
    for j, smp in enumerate(truth.columns):
        w, _ = nnls(A, Y[:, j]); w = w / w.sum() if w.sum() > 0 else w
        est = pd.Series(w, index=types).groupby(lambda t: GROUP[t]).sum().reindex(groups)
        E.append(est.values)
        for gname, v, tv in zip(groups, est.values, truth[smp].values):
            est_rows.append(dict(strategy=name, sample=smp, group=gname, estimate=float(v), flow=float(tv)))
    E = np.array(E).T                                       # groups x samples
    T = truth.values
    err = E - T
    met_rows.append(dict(strategy=name, rmse=float(np.sqrt((err ** 2).mean())), r=float(np.corrcoef(E.ravel(), T.ravel())[0, 1]),
                         rmse_per_sample_min=float(np.sqrt((err ** 2).mean(0)).min()),
                         rmse_per_sample_max=float(np.sqrt((err ** 2).mean(0)).max()),
                         **{f"bias_{g}": float(err[i].mean()) for i, g in enumerate(groups)}))
    log(f"{name:15} RMSE {met_rows[-1]['rmse']:.3f}  r {met_rows[-1]['r']:+.2f}  " +
        " ".join(f"{g}:{err[i].mean():+.3f}" for i, g in enumerate(groups)))
pd.DataFrame(est_rows).to_csv(HERE / "real_bulk_estimates.csv", index=False)
pd.DataFrame(met_rows).to_csv(HERE / "real_bulk_metrics.csv", index=False)
med = pd.Series({t: np.median(np.asarray(X[lab == t].sum(1)).ravel()) for t in types})
log("median reference row sum per type: " + ", ".join(f"{t} {v:.0f}" for t, v in med.items()))
log("done")
