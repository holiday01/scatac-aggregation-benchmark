"""
02: The six aggregation strategies of the paper, applied unchanged to the adult
    non-tumor tissues of the human chromatin accessibility atlas (GSE184462).

The aggregation functions and the permutation estimator are imported from
tools/benchmark_aggregation.py, so every strategy is the same code that produced the
tumor results. Cell types are the atlas's own labels; types with fewer than
MIN_CELLS cells in a tissue are dropped.

Read-outs per (tissue, model, strategy), as in the paper:
  1. top share: largest share of per-gene maxima held by one cell type, and whether
     that type is the one with the highest median per-cell coverage;
  2. Spearman rho across types between median coverage (n_fragment) and share of maxima;
  3. lineage-marker concordance on broad lineages shared by most tissues. A marker is
     correct when its maximum falls in any atlas type that belongs to its lineage.
     The atlas labels were assigned from gene activity, so this panel is NOT
     independent of the labels and is reported as a secondary read-out only.
  4. leave-one-sample-out range of the top share (tissues with >= 3 samples).

Run:  python3 02_benchmark_catlas.py                 # every tissue with built matrices
      python3 02_benchmark_catlas.py liver lung      # selected tissues
Out:  catlas_benchmark.csv, catlas_depth_by_type.csv, catlas_markers.csv, catlas_loso.csv
"""
import sys, re, time
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
import benchmark_aggregation as B  # noqa: E402

DATA = Path("/mnt/10t/scrna_atac/data/processed/catlas_adult")
MIN_CELLS = 100
SEED = 42

# broad lineage of an atlas cell-type name (first match wins); unmatched types stay in
# the analysis but no marker targets them
LINEAGE = [
    (r"Macrophage", "Macrophage"),
    (r"T Lymphocyte|Natural Killer", "T_NK"),
    (r"Lymphatic Endothelial", None),
    (r"Endothelial|Endocardial", "Endothelial"),
    (r"Fibroblast", "Fibroblast"),
    (r"Smooth Muscle", "Smooth_muscle"),
]
MARKERS = {
    "Macrophage":    ["MS4A7", "FOLR2", "CD163", "C1QA", "C1QB", "MPEG1"],
    "T_NK":          ["CD2", "CD247", "SKAP1", "THEMIS", "CD96", "ITK"],
    "Endothelial":   ["CLDN5", "FLT1", "ERG", "ENG", "TEK", "PLVAP"],
    "Fibroblast":    ["COL3A1", "PDGFRA", "FBLN1", "COL6A3", "COL5A1", "DCN"],
    "Smooth_muscle": ["MYH11", "CNN1", "LMOD1", "DES", "ACTG2", "SYNPO2"],
}


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def lineage_of(name):
    for pat, lin in LINEAGE:
        if re.search(pat, name): return lin
    return None


def top_share(P, expressed):
    arg = P[:, expressed].argmax(0)
    share = np.bincount(arg, minlength=P.shape[0]) / arg.size
    return share


def strategies(X, lab, types, rng, counts):
    out = [("raw_mean", lambda: B.m_raw(X, lab, types)),
           ("zscale_mean", lambda: B.m_zscale(X, lab, types)),
           ("lognorm_mean", lambda: B.m_lognorm(X, lab, types)),
           ("pseudobulk_cpm", lambda: B.m_pseudobulk_cpm(X, lab, types)),
           ("cpm_mean", lambda: B.m_cpm(X, lab, types))]
    if counts:  # binomial thinning is defined for integer counts only
        out.append(("depth_matched", lambda: B.m_depth_matched(X, lab, types, rng)[0]))
    return out


def main():
    organs = sys.argv[1:] or sorted({p.name.split("_geneactivity_")[0] for p in DATA.glob("*_geneactivity_proximal.h5ad")})
    rng = np.random.default_rng(SEED)
    rows, drows, mrows, lrows = [], [], [], []
    for organ in organs:
        for model in ("proximal", "enhancer"):
            path = DATA / f"{organ}_geneactivity_{model}.h5ad"
            if not path.exists(): continue
            a = ad.read_h5ad(path)
            n_by = a.obs["cell_type"].value_counts()
            a = a[a.obs["cell_type"].isin(n_by.index[n_by >= MIN_CELLS])].copy()
            lab = a.obs["cell_type"].astype(str).values
            types = sorted(pd.unique(lab))
            if len(types) < 3:
                log(f"{organ}/{model}: {len(types)} types with >= {MIN_CELLS} cells, skipped"); continue
            X = a.X.tocsr().astype(np.float32)
            genes = a.var_names.astype(str).values; gidx = {g: i for i, g in enumerate(genes)}
            counts = bool(np.allclose(X.data[:100000], np.rint(X.data[:100000])))
            dep = a.obs["n_fragment"].values.astype(float)
            sample = a.obs["sample"].astype(str).values
            med = np.array([np.median(dep[lab == t]) for t in types])
            spread = float(med.max() / med.min())
            deepest = types[int(med.argmax())]
            expressed = np.asarray(X.sum(0)).ravel() > 0
            log(f"{organ}/{model}: {X.shape[0]} cells, {len(types)} types, {len(set(sample))} samples, "
                f"coverage spread {spread:.2f}x, deepest {deepest}")
            for t, m_ in zip(types, med):
                drows.append(dict(organ=organ, model=model, cell_type=t, lineage=lineage_of(t),
                                  n_cells=int((lab == t).sum()), median_n_fragment=float(m_)))
            lin = {t: lineage_of(t) for t in types}
            present = {l for l in lin.values() if l}
            panel = [(g, l) for l, gs in MARKERS.items() if l in present for g in gs if g in gidx]

            for name, fn in strategies(X, lab, types, rng, counts):
                P = fn()
                share = top_share(P, expressed); top = int(share.argmax())
                rho, _ = spearmanr(med, share)
                marg = np.array([lin[types[int(P[:, gidx[g]].argmax())]] or "other" for g, _ in panel])
                mtgt = np.array([l for _, l in panel])
                obs = int((marg == mtgt).sum()) if len(panel) else 0
                pp, nmean, _ = B.perm_p(marg, mtgt, obs, rng) if len(panel) else (np.nan, np.nan, np.nan)
                rows.append(dict(organ=organ, model=model, method=name, n_cells=X.shape[0], n_types=len(types),
                                 n_samples=len(set(sample)), coverage_spread=spread,
                                 top_share=float(share[top]), uniform_share=1.0 / len(types),
                                 top_type=types[top], deepest_type=deepest, top_is_deepest=bool(types[top] == deepest),
                                 rho_coverage_share=float(rho),
                                 markers_correct=obs, markers_total=len(panel), perm_p=pp, null_mean=nmean))
                for (g, l), ag in zip(panel, marg):
                    mrows.append(dict(organ=organ, model=model, method=name, gene=g, target=l, argmax_lineage=ag,
                                      correct=int(ag == l)))
                log(f"  {name:15} top {types[top][:32]:>32} {share[top]:5.1%} (1/k {1/len(types):.1%})  "
                    f"rho={rho:+.2f}  markers {obs}/{len(panel)}")

            # leave-one-sample-out, the three strategies shown in Fig 1D
            ss = sorted(set(sample))
            if len(ss) >= 3:
                for s in ss:
                    k = sample != s
                    Xk, labk = X[k], lab[k]
                    tk = [t for t in types if (labk == t).sum() >= 20]
                    ek = np.asarray(Xk.sum(0)).ravel() > 0
                    for name, fn in (("raw_mean", B.m_raw), ("lognorm_mean", B.m_lognorm), ("cpm_mean", B.m_cpm)):
                        sh = top_share(fn(Xk, labk, tk), ek)
                        lrows.append(dict(organ=organ, model=model, method=name, omitted=s,
                                          n_cells=int(k.sum()), top_share=float(sh.max()), top_type=tk[int(sh.argmax())]))
            del a, X
        tag = "" if not sys.argv[1:] else "_partial"
        pd.DataFrame(rows).to_csv(HERE / f"catlas_benchmark{tag}.csv", index=False)
        pd.DataFrame(drows).to_csv(HERE / f"catlas_depth_by_type{tag}.csv", index=False)
        pd.DataFrame(mrows).to_csv(HERE / f"catlas_markers{tag}.csv", index=False)
        pd.DataFrame(lrows).to_csv(HERE / f"catlas_loso{tag}.csv", index=False)

    R = pd.DataFrame(rows)
    print("\n=== median top share across tissues, by model and strategy ===")
    for model in ("proximal", "enhancer"):
        for name in ["raw_mean", "zscale_mean", "lognorm_mean", "pseudobulk_cpm", "cpm_mean", "depth_matched"]:
            s = R[(R.model == model) & (R.method == name)]
            if len(s) == 0: continue
            print(f"  {model:9} {name:15} {s.top_share.median():5.1%} ({s.top_share.min():.0%}-{s.top_share.max():.0%})  "
                  f"top type is deepest {int(s.top_is_deepest.sum())}/{len(s)}  rho median {s.rho_coverage_share.median():+.2f}  "
                  f"markers {int(s.markers_correct.sum())}/{int(s.markers_total.sum())}")


if __name__ == "__main__":
    main()
