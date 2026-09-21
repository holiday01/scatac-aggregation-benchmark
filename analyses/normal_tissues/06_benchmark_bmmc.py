"""
06: Aggregation benchmark on a THIRD-PARTY gene-activity matrix with a multi-donor RNA reference:
    the NeurIPS 2021 open-problems bone-marrow multiome (GSE194122; 69,249 cells, 10 healthy
    donors, 13 batches, 4 sites, 22 expert-annotated types).

Gene activity is the matrix deposited by the authors (obsm['ATAC_gene_activity'], Signac
GeneActivity), not one built here. It is deposited log-normalised, log1p(count * s / total);
the integer counts are recovered exactly by dividing expm1(value) by each cell's smallest
non-zero expm1(value) (every entry lands within 1e-4 of an integer; asserted below).
Per-cell coverage is obs['ATAC_atac_fragments'] (never a matrix row sum).

Read-outs, per strategy (functions imported from tools/benchmark_aggregation.py):
  top share / modal-argmax type / rho(median coverage, share of maxima)
  agreement of the per-gene maximum with the RNA reference (CP10K mean of GEX counts),
    over all shared genes and over RNA-informative genes (top type >= 2x the second)
  the same, leaving one donor out at a time (raw, log-normalised, CP10K, as deposited)

Out: bmmc_benchmark.csv, bmmc_depth_by_type.csv (with the RNA reference's own share), bmmc_lodo.csv
"""
import sys, time, warnings
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "tools"))
from benchmark_aggregation import type_mean, m_raw, m_zscale, m_pseudobulk_cpm, m_depth_matched, m_cpm, m_lognorm  # noqa

H5AD = Path("/mnt/10t/scrna_atac/data/raw/multiome_bmmc_neurips/bmmc_multiome.h5ad")
MIN_CELLS, SEED = 100, 42


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


a = ad.read_h5ad(H5AD)
obs = a.obs
lab_all = obs["cell_type"].astype(str).values
vc = pd.Series(lab_all).value_counts()
types = sorted(vc.index[vc >= MIN_CELLS])
dep = obs["ATAC_atac_fragments"].values.astype(float)
donor = obs["DonorID"].astype(str).values
log(f"{a.shape[0]} cells, {len(types)} types with >= {MIN_CELLS} cells, {len(set(donor))} donors, {obs['batch'].nunique()} batches")

# ── deposited gene activity -> exact integer counts ──────────────────────────
D = sp.csr_matrix(a.obsm["ATAC_gene_activity"]).astype(np.float64)       # as deposited (log-normalised)
genes = np.asarray(a.uns["ATAC_gene_activity_var_names"]).astype(str)
E = D.copy(); E.data = np.expm1(E.data)
n = np.diff(E.indptr); assert (n > 0).all()
unit = np.minimum.reduceat(E.data, E.indptr[:-1])
C = E.copy(); C.data = C.data / np.repeat(unit, n)
dev = np.abs(C.data - np.rint(C.data)).max()
assert dev < 1e-3, dev
C.data = np.rint(C.data); X = C.astype(np.float32).tocsr(); D = D.astype(np.float32).tocsr()
log(f"recovered counts: max deviation from integer {dev:.1e}; row-sum median {np.median(np.asarray(X.sum(1)).ravel()):.0f}; "
    f"r(row sum, fragments) = {np.corrcoef(np.asarray(X.sum(1)).ravel(), dep)[0, 1]:.3f}")

# ── RNA reference: CP10K mean of GEX counts ──────────────────────────────────
gex = (a.var["feature_types"] == "GEX").values
Xr = sp.csr_matrix(a.layers["counts"][:, gex]).astype(np.float32)
rna_genes = a.var_names[gex].astype(str).values
del a


def rna_reference(keep):
    R = m_cpm(Xr[keep], lab_all[keep], types)
    det = np.asarray((Xr[keep] > 0).mean(0)).ravel()
    srt = np.sort(R, 0)
    return R, det, (srt[-1] >= 2 * np.maximum(srt[-2], 1e-9)) & (srt[-1] > 0)


def readout(P, expressed, med, ai, r_arg, inf):
    arg = P[:, expressed].argmax(0)
    share = np.bincount(arg, minlength=len(types)) / arg.size
    a_arg = P[:, ai].argmax(0)
    return share, float(spearmanr(med, share)[0]), float((a_arg == r_arg).mean()), float((a_arg[inf] == r_arg[inf]).mean())


def shared_genes(expressed, det):
    rmap = {g: i for i, g in enumerate(rna_genes)}
    pairs = [(i, rmap[g]) for i, g in enumerate(genes) if g in rmap and expressed[i] and det[rmap[g]] >= 0.01]
    return np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])


rng = np.random.default_rng(SEED)
keep = np.isin(lab_all, types)
lab = lab_all[keep]
STRATS = [("raw_mean", lambda M, Dm, l: m_raw(M, l, types)),
          ("zscale_mean", lambda M, Dm, l: m_zscale(M, l, types)),
          ("lognorm_mean", lambda M, Dm, l: m_lognorm(M, l, types)),
          ("deposited_mean", lambda M, Dm, l: type_mean(Dm, l, types)),     # what averaging the deposited matrix gives
          ("pseudobulk_cpm", lambda M, Dm, l: m_pseudobulk_cpm(M, l, types)),
          ("cpm_mean", lambda M, Dm, l: m_cpm(M, l, types)),
          ("depth_matched", lambda M, Dm, l: m_depth_matched(M, l, types, rng)[0])]

# ── all donors ───────────────────────────────────────────────────────────────
Xk, Dk = X[keep], D[keep]
med = np.array([np.median(dep[keep][lab == t]) for t in types])
expressed = np.asarray(Xk.sum(0)).ravel() > 0
R, det, inf_all = rna_reference(keep)
ai, ri = shared_genes(expressed, det)
r_arg, inf = R[:, ri].argmax(0), inf_all[ri]
rna_share = np.bincount(r_arg, minlength=len(types)) / len(r_arg)
log(f"coverage spread {med.max() / med.min():.2f}x (min {types[int(med.argmin())]} {med.min():.0f}, max {types[int(med.argmax())]} {med.max():.0f}); "
    f"shared genes {len(ai)}, RNA-informative {int(inf.sum())}; RNA reference top share {rna_share.max():.1%} ({types[int(rna_share.argmax())]})")
pd.DataFrame(dict(cell_type=types, n_cells=[int((lab == t).sum()) for t in types], median_fragments=med,
                  rna_reference_share=rna_share)).to_csv(HERE / "bmmc_depth_by_type.csv", index=False)
rows = []
for name, fn in STRATS:
    t0 = time.time()
    share, rho, ag, agi = readout(fn(Xk, Dk, lab), expressed, med, ai, r_arg, inf)
    top = int(share.argmax())
    rows.append(dict(method=name, n_cells=int(keep.sum()), n_types=len(types), coverage_spread=float(med.max() / med.min()),
                     top_share=float(share[top]), top_type=types[top], deepest_type=types[int(med.argmax())],
                     rho_coverage_share=rho, share_l1_to_rna=float(np.abs(share - rna_share).sum()),
                     argmax_agree_all=ag, argmax_agree_informative=agi, n_shared_genes=len(ai), n_informative=int(inf.sum())))
    log(f"  {name:15} top {types[top]:>16} {share[top]:5.1%}  rho={rho:+.2f}  RNA agree {ag:.3f} (inf {agi:.3f})  "
        f"L1 to RNA share {np.abs(share - rna_share).sum():.2f}  ({time.time() - t0:.0f}s)")
pd.DataFrame(rows).to_csv(HERE / "bmmc_benchmark.csv", index=False)

# ── leave one donor out ──────────────────────────────────────────────────────
lrows = []
for d in sorted(set(donor)):
    k = keep & (donor != d)
    labd = lab_all[k]
    if min((labd == t).sum() for t in types) < 20:
        log(f"  donor {d}: a type falls below 20 cells, fold skipped"); continue
    Xd, Dd = X[k], D[k]
    medd = np.array([np.median(dep[k][labd == t]) for t in types])
    ex = np.asarray(Xd.sum(0)).ravel() > 0
    Rd, detd, infd = rna_reference(k)
    aid, rid = shared_genes(ex, detd)
    for name, fn in STRATS:
        if name in ("zscale_mean", "pseudobulk_cpm", "depth_matched"): continue
        share, rho, ag, agi = readout(fn(Xd, Dd, labd), ex, medd, aid, Rd[:, rid].argmax(0), infd[rid])
        lrows.append(dict(omitted_donor=d, method=name, n_cells=int(k.sum()), top_share=float(share.max()),
                          top_type=types[int(share.argmax())], rho_coverage_share=rho,
                          argmax_agree_all=ag, argmax_agree_informative=agi))
    log(f"  donor {d} omitted ({int(k.sum())} cells): " + "  ".join(
        f"{r['method']} {r['top_share']:.1%}/{r['argmax_agree_informative']:.3f}" for r in lrows[-4:]))
pd.DataFrame(lrows).to_csv(HERE / "bmmc_lodo.csv", index=False)
L = pd.DataFrame(lrows)
print("\n=== leave-one-donor-out ranges ===")
print(L.groupby("method").agg(top_min=("top_share", "min"), top_max=("top_share", "max"),
                              agree_inf_min=("argmax_agree_informative", "min"),
                              agree_inf_max=("argmax_agree_informative", "max")).round(3).to_string())
