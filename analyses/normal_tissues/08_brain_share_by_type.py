"""
08: Share of per-gene maxima held by every cell type in the brain multiome, per aggregation strategy.

brain3k_benchmark.csv (05) stores only the largest share and its type. The manuscript also quotes the
share held by excitatory neurons when that type is not the largest (CP10K mean, weighted model), so
this script writes the full per-type distribution. Same cells, labels, matrices and aggregation
functions as 05; the coverage-matched control is omitted because its thinning draws depend on the
random-number sequence of 05 (its top share is in brain3k_benchmark.csv).

Out: brain3k_share_by_type.csv  (model, method, cell_type, share_of_maxima, n_genes_expressed)
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_aggregation import m_raw, m_zscale, m_pseudobulk_cpm, m_cpm, m_lognorm, m_archr  # noqa

HERE = Path(__file__).resolve().parent
PROC = Path("/mnt/10t/scrna_atac/data/processed/multiome_brain3k")
MODELS = [("proximal", PROC / "brain3k_geneactivity_proximal.h5ad"),
          ("enhancer", PROC / "brain3k_geneactivity_enhancer.h5ad")]
METHODS = [("raw_mean", m_raw), ("zscale_mean", m_zscale), ("lognorm_mean", m_lognorm),
           ("pseudobulk_cpm", m_pseudobulk_cpm), ("cpm_mean", m_cpm), ("archr_style", m_archr)]

lab_df = pd.read_csv(HERE / "brain3k_rna_labels.csv", index_col=0)
lab_df = lab_df[lab_df.cell_type != "dropped"]
types = sorted(lab_df.cell_type.unique())
bench = pd.read_csv(HERE / "brain3k_benchmark.csv").set_index(["model", "method"])

rows = []
for model, path in MODELS:
    a = ad.read_h5ad(path)
    a = a[lab_df.index.intersection(a.obs_names)].copy()
    lab = lab_df.loc[a.obs_names, "cell_type"].values.astype(str)
    X = sp.csr_matrix(a.X).astype(np.float32).tocsr()
    expressed = np.asarray(X.sum(0)).ravel() > 0
    for name, fn in METHODS:
        res = fn(X, lab, types); P = res[0] if isinstance(res, tuple) else res
        arg = P[:, expressed].argmax(0)
        share = np.bincount(arg, minlength=len(types)) / arg.size
        # consistency with 05: same largest share and type
        b = bench.loc[(model, name)]
        assert abs(share.max() - b.top_share) < 1e-12 and types[int(share.argmax())] == b.top_type, (model, name)
        for t, s in zip(types, share):
            rows.append(dict(model=model, method=name, cell_type=t, share_of_maxima=float(s),
                             n_genes_expressed=int(expressed.sum())))
        print(model, name, " ".join(f"{t}={s:.4f}" for t, s in zip(types, share)), flush=True)

pd.DataFrame(rows).to_csv(HERE / "brain3k_share_by_type.csv", index=False)
print("wrote brain3k_share_by_type.csv")
