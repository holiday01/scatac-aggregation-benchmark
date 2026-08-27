"""
Script 91: Does a CORRECT chromatin prior change the result?
============================================================
Script 90 rebuilt the HCC scATAC gene activity from fragments with real
coordinates, in two flavours:
    proximal  = gene body + 2 kb upstream   (== Signac GeneActivity, the model
                Craig et al. actually used for GSE227265) — correctly named
    enhancer  = ArchR-style distance-weighted score over 5 kb tiles within
                100 kb (exp decay, tau=5 kb) — sees distal enhancers

The published prior (`atac_hcc_geneactivity.h5ad`) has MISASSIGNED gene names
(Script 04 guessed the row order), so the ATAC arm of HCC and HNSC — 6 of the
paper's 9 combinations — was selected with a chromatin prior carrying little
valid lineage signal.

This script holds EVERYTHING ELSE constant (Script 76's formulas, N_POOL=100,
N_EQUAL=25; Script 83's paired bootstrap C-index test) and swaps only the prior:

    old_guessed | new_proximal | new_enhancer

Two questions, finally answerable:
  Q1 (naming)  old_guessed vs new_proximal — was the null an artefact of the bug?
  Q2 (enhancer) new_proximal vs new_enhancer — does an enhancer-aware prior help?

Outputs:
  results/validation/corrected_prior_markers.csv
  results/validation/corrected_prior_headtohead.csv
"""
import sys
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from scipy import stats
import os
from pathlib import Path
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

BASE    = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
OUT_DIR = BASE / "results/validation"
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_POOL, N_EQUAL = 100, 25          # identical to Script 76
N_BOOT, RNG_SEED = 2000, 42        # identical to Script 83

PRIORS = {
    'old_guessed':  BASE / "data/processed/integrated/atac_hcc_geneactivity.h5ad",
    'new_proximal': BASE / "data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad",
    'new_enhancer': BASE / "data/processed/integrated/atac_hcc_geneactivity_enhancer.h5ad",
}

# HCC and HNSC both use the HCC scATAC reference (Script 76). COAD has its own
# (separately broken: z-scored) and is out of scope here.
CONFIGS = {
    'HCC':  {'scrna_ref': BASE / 'data/processed/integrated/pseudobulk_reference_v2.csv',
             'expr': BASE / 'data/raw/TCGA_LIHC/TCGA_LIHC_expression.gz',
             'clin': BASE / 'data/raw/TCGA_LIHC/TCGA_LIHC_clinical.tsv',
             'pairs': [('Macrophage', 'Macrophage'), ('B_cell', 'B_cell'), ('NK', 'NK_cytotoxic_T')]},
    'HNSC': {'scrna_ref': BASE / 'data/processed/HNSC/pseudobulk_reference_hnsc_v3.csv',
             'expr': BASE / 'data/raw/TCGA_HNSC/TCGA_HNSC_expression.gz',
             'clin': BASE / 'data/raw/TCGA_HNSC/TCGA_HNSC_clinical.tsv',
             'pairs': [('Macrophage', 'Macrophage'), ('B_Plasma', 'B_cell'), ('CD8_T', 'NK_cytotoxic_T')]},
}


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


NORMALISE = '--raw' not in sys.argv   # --raw reproduces the superseded numbers


def atac_celltype_mean(h5):
    """Mean gene activity per cell type, depth-normalised per cell first.

    Without the normalisation the group mean is a sum over raw counts, so the
    deepest cell type wins the argmax for genes unrelated to it -- on these priors
    the B-cell compartment takes 66-76 % of all genome-wide argmaxima raw, against
    ~15 % normalised. `atac_spec` below is a ratio to `am.max(1)`, so it degenerates
    to a constant for whichever type dominates, and the rank product then carries no
    chromatin information at all for that lineage.
    """
    a = ad.read_h5ad(h5)
    X = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    X = X.astype(np.float64)
    if NORMALISE:
        tot = X.sum(1, keepdims=True)
        tot[tot == 0] = 1.0
        X = X / tot * 1e4
    ct = a.obs['cell_type'].astype(str).values
    return pd.DataFrame({c: X[ct == c].mean(0) for c in pd.unique(ct)}, index=a.var_names)


def derive_markers(ref, atac_mean, rna_ct, atac_ct):
    """Script 76 logic, verbatim: scRNA specificity -> pool -> rank product."""
    ref_c = ref[[c for c in ref.columns]]
    others = [c for c in ref_c.columns if c != rna_ct]
    spec = ref_c[rna_ct] / (ref_c[others].mean(1) + 0.01)      # scRNA specificity
    pool = spec.nlargest(N_POOL).index

    scrna_top = list(spec.nlargest(N_EQUAL).index)             # scRNA-only arm

    common = [g for g in pool if g in atac_mean.index]
    if atac_ct not in atac_mean.columns or len(common) < N_EQUAL:
        return None, scrna_top
    am = atac_mean.loc[common]
    atac_spec = am[atac_ct] / (am.max(1) + 1e-6)               # ATAC specificity
    r_rna  = spec.loc[common].rank(ascending=False)
    r_atac = atac_spec.rank(ascending=False)
    rank_prod = r_rna * r_atac
    atac_top = list(rank_prod.nsmallest(N_EQUAL).index)        # ATAC-informed arm
    return atac_top, scrna_top


def parse_os(clin):
    ev = clin.get('vital_status', pd.Series('', index=clin.index)).map(
        lambda x: 1 if str(x).upper() in ('DEAD', 'DECEASED', '1') else
                  (0 if str(x).upper() in ('ALIVE', 'LIVING', '0') else np.nan))
    dtd  = pd.to_numeric(clin.get('days_to_death', pd.Series(dtype=float, index=clin.index)), errors='coerce')
    dtlf = pd.to_numeric(clin.get('days_to_last_followup',
           clin.get('days_to_last_follow_up', pd.Series(dtype=float, index=clin.index))), errors='coerce')
    return pd.to_numeric(ev, errors='coerce'), pd.to_numeric(dtd.fillna(dtlf), errors='coerce')


def c_index(t, score, ev):
    from lifelines.utils import concordance_index
    return concordance_index(t, -score, ev)


def norm_bc(b):
    p = str(b).split('-')
    return '-'.join(p[:4]) if len(p) >= 4 else str(b)[:15]


def load_cohort(cfg):
    expr = pd.read_csv(cfg['expr'], sep='\t', index_col=0)
    clin = pd.read_csv(cfg['clin'], sep='\t', index_col=0)
    expr.columns = [norm_bc(c) for c in expr.columns]
    clin.index   = [norm_bc(i) for i in clin.index]
    ev, t = parse_os(clin)
    common = [c for c in expr.columns if c in clin.index]
    return expr[common], ev.reindex(common), t.reindex(common)


def mean_score(expr, genes):
    avail = [g for g in genes if g in expr.index]
    if len(avail) < 3:
        return None
    return np.log1p(expr.loc[avail]).mean(axis=0)


# ── run ──────────────────────────────────────────────────────────────────────
atac_means = {}
for tag, h5 in PRIORS.items():
    if not h5.exists():
        log(f"MISSING prior {tag}: {h5} — run Script 90 first"); continue
    atac_means[tag] = atac_celltype_mean(h5)
    log(f"loaded prior {tag}: {atac_means[tag].shape[0]} genes x {atac_means[tag].shape[1]} cell types")

marker_rows, res_rows = [], []

for cancer, cfg in CONFIGS.items():
    log(f"=== {cancer} ===")
    ref = pd.read_csv(cfg['scrna_ref'], index_col=0)
    expr, ev, t = load_cohort(cfg)

    for rna_ct, atac_ct in cfg['pairs']:
        if rna_ct not in ref.columns:
            log(f"  {rna_ct}: not in scRNA ref — skipped"); continue

        for tag, am in atac_means.items():
            atac_top, scrna_top = derive_markers(ref, am, rna_ct, atac_ct)
            if atac_top is None:
                log(f"  {cancer}/{rna_ct}/{tag}: prior lacks cell type — skipped"); continue

            for g in atac_top:
                marker_rows.append(dict(cancer=cancer, cell_type=rna_ct, prior=tag,
                                        arm='ATAC-informed', gene=g))

            sa = mean_score(expr, atac_top)
            ss = mean_score(expr, scrna_top)
            if sa is None or ss is None:
                continue
            d = pd.DataFrame({'a': sa, 's': ss, 'ev': ev, 't': t}).dropna()
            d = d[d['t'] > 0]
            if len(d) < 30 or d['ev'].sum() < 10:
                continue
            tt, ee = d['t'].values, d['ev'].values
            a, s = d['a'].values, d['s'].values

            Ca, Cs = c_index(tt, a, ee), c_index(tt, s, ee)
            dAbs = abs(Ca - 0.5) - abs(Cs - 0.5)

            rng = np.random.default_rng(RNG_SEED)
            bd = []
            n = len(d)
            for _ in range(N_BOOT):
                i = rng.integers(0, n, n)
                if ee[i].sum() < 5:
                    continue
                ca = c_index(tt[i], a[i], ee[i]); cs = c_index(tt[i], s[i], ee[i])
                bd.append(abs(ca - 0.5) - abs(cs - 0.5))
            bd = np.array(bd)
            lo, hi = np.percentile(bd, [2.5, 97.5])
            p2 = float(2 * min((bd <= 0).mean(), (bd >= 0).mean()))

            n_overlap = len(set(atac_top) & set(scrna_top))
            res_rows.append(dict(cancer=cancer, cell_type=rna_ct, prior=tag,
                                 n=n, events=int(ee.sum()),
                                 C_atac=round(Ca, 4), C_scrna=round(Cs, 4),
                                 dAbsC=round(dAbs, 4),
                                 ci_lo=round(lo, 4), ci_hi=round(hi, 4),
                                 p_twosided=round(min(p2, 1.0), 4),
                                 overlap_with_scrna_arm=n_overlap))
            log(f"  {rna_ct:11} {tag:13} C_atac={Ca:.3f} C_scrna={Cs:.3f} "
                f"dAbsC={dAbs:+.4f} p={min(p2,1.0):.3f} overlap={n_overlap}/25")

mk = pd.DataFrame(marker_rows); mk.to_csv(OUT_DIR / "corrected_prior_markers.csv", index=False)
rs = pd.DataFrame(res_rows);    rs.to_csv(OUT_DIR / "corrected_prior_headtohead.csv", index=False)

print("\n" + "=" * 70)
print("CROSS-COMBINATION PAIRED TEST, BY PRIOR  (H1: ATAC-informed better)")
print("=" * 70)
for tag in PRIORS:
    sub = rs[rs.prior == tag]
    if len(sub) < 2:
        continue
    d = sub.dAbsC.values
    try:
        w = stats.wilcoxon(d, alternative='greater', zero_method='wilcox').pvalue
    except Exception:
        w = float('nan')
    sg = stats.binomtest(int((d > 0).sum()), len(d), 0.5, alternative='greater').pvalue
    print(f"  {tag:13}  n={len(d)}  ATAC better in {int((d>0).sum())}/{len(d)}  "
          f"mean dAbsC={d.mean():+.4f}  Wilcoxon p={w:.4f}  sign p={sg:.4f}")
print()
log("Script 91 done.")
