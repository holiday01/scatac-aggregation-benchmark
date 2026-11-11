"""Annotation task over repeated stratified splits.

The single 50/50 split used in 03_benchmark_multiome.py gives one number per
strategy with no indication of split-to-split variability, and Pearson
correlation is invariant to a positive rescaling of a reference row, so the raw
mean and pooled CPM are identical by construction under it. This script repeats
the task over five stratified splits and reports both Pearson and rank
(Spearman) accuracy as mean and range.

Out: multiome_annotation_repeats.csv (per split) and a printed summary.
"""
import sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
from scipy.stats import rankdata
warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_aggregation import (m_raw, m_zscale, m_pseudobulk_cpm, m_depth_matched,
                                   m_cpm, m_lognorm, m_archr, percell_scale)  # noqa

HERE = Path(__file__).resolve().parent
PROC = Path("/mnt/10t/scrna_atac/data/processed/multiome_pbmc10k")
MODELS = [("proximal", PROC / "pbmc10k_geneactivity_proximal.h5ad"),
          ("enhancer", PROC / "pbmc10k_geneactivity_enhancer.h5ad")]
SPLITS = 5
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)

lab_df = pd.read_csv(HERE / "rna_labels.csv", index_col=0)
lab_df = lab_df[lab_df.cell_type != "dropped"]
types = sorted(lab_df.cell_type.unique())

def corr_assign(Q, P):
    Q = Q - Q.mean(1, keepdims=True); P = P - P.mean(1, keepdims=True)
    Q = Q / (np.linalg.norm(Q, axis=1, keepdims=True) + 1e-12)
    P = P / (np.linalg.norm(P, axis=1, keepdims=True) + 1e-12)
    return (Q @ P.T).argmax(1)

rows = []
for model, path in MODELS:
    a = ad.read_h5ad(path)
    keep = a.obs_names.isin(lab_df.index)
    X = (a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)).tocsr().astype(np.float32)[keep]
    lab = lab_df.cell_type.reindex(a.obs_names[keep]).values.astype(str)
    log(f"{model}: {X.shape[0]} cells, {X.shape[1]} genes")
    for seed in range(SPLITS):
        rng = np.random.default_rng(100 + seed)
        train = np.zeros(X.shape[0], bool)
        for t in types:
            idx = np.where(lab == t)[0]; rng.shuffle(idx); train[idx[: len(idx) // 2]] = True
        test = ~train
        Xn = percell_scale(X[test]).tocsr(); Xn.data = np.log1p(Xn.data)
        genes = np.asarray(X[train].sum(0)).ravel() > 0
        Q_log = Xn[:, genes].toarray().astype(np.float64)
        Q_rank = np.vstack([rankdata(r) for r in Q_log]).astype(np.float64)
        ylab = lab[test]; Xtr, labtr = X[train], lab[train]
        for name, fn in [("raw_mean", m_raw), ("zscale_mean", m_zscale), ("lognorm_mean", m_lognorm),
                         ("pseudobulk_cpm", m_pseudobulk_cpm), ("cpm_mean", m_cpm), ("archr_style", m_archr),
                         ("depth_matched", lambda M, l, T: m_depth_matched(M, l, T, rng))]:
            P = fn(Xtr, labtr, types)
            if isinstance(P, tuple): P = P[0]
            Pm = P[:, genes].astype(np.float64)
            pred_p = corr_assign(Q_log.copy(), Pm.copy())
            Pr = np.vstack([rankdata(r) for r in Pm]).astype(np.float64)
            pred_s = corr_assign(Q_rank.copy(), Pr)
            r = dict(model=model, split=seed, method=name, n_test=int(test.sum()))
            for tag, pred in [("pearson", pred_p), ("spearman", pred_s)]:
                y = np.array(types)[pred]
                r[f"acc_{tag}"] = float((y == ylab).mean())
                r[f"balacc_{tag}"] = float(np.mean([(y[ylab == t] == t).mean() for t in types]))
            rows.append(r)
        log(f"  split {seed}: done")
    del a, X

D = pd.DataFrame(rows)
D.to_csv(HERE / "multiome_annotation_repeats.csv", index=False)
print("\n=== annotation over %d stratified splits (mean [min-max]) ===" % SPLITS)
for model in ("proximal", "enhancer"):
    print(f"-- {model}")
    for m in D.method.unique():
        s = D[(D.model == model) & (D.method == m)]
        print(f"   {m:15} " + "  ".join(
            f"{tag} {s[c].mean():.3f} [{s[c].min():.3f}-{s[c].max():.3f}]"
            for tag, c in [("Pearson acc", "acc_pearson"), ("rank acc", "acc_spearman"),
                           ("rank balacc", "balacc_spearman")]))
print(f"\nwrote {HERE/'multiome_annotation_repeats.csv'}")
