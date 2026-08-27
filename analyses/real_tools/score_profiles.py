"""
Shared read-outs for the real-tool runs, copied verbatim in logic from
tools/benchmark_aggregation.py so that the numbers are comparable:
  - top argmax share over EXPRESSED genes (column sum > 0 in the counts matrix)
  - held-out marker concordance on PANEL['HCC'] with the 10,000-shuffle permutation null
Profiles are (types x genes) arrays; `genes` must be aligned to the counts matrix columns.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tools"))
from benchmark_aggregation import PANEL, true_depth, BASE, perm_p, N_PERM, SEED  # noqa

SCRATCH = Path("/tmp/claude-1000/-mnt-10t-scrna-atac/b849b398-690b-4bff-b242-236b530a2a05/scratchpad")
OUTDIR = Path(__file__).resolve().parent
RELABEL = {"B_cell": "Endothelial_stromal", "DC": "B_cell"}


def relabel(lab):
    lab = np.asarray(lab).astype(str)
    return np.where(lab == "B_cell", "Endothelial_stromal", np.where(lab == "DC", "B_cell", lab))


def score(P, types, genes, expressed, cohort="HCC", seed=SEED):
    """P: (n_types x n_genes) profile. Returns dict of read-outs + per-marker table."""
    rng = np.random.default_rng(seed)
    P = np.asarray(P, dtype=np.float64)
    types = list(types); genes = np.asarray(genes).astype(str)
    arg = P[:, expressed].argmax(0)
    share = np.bincount(arg, minlength=len(types)) / arg.size
    top = int(share.argmax())
    gidx = {g: i for i, g in enumerate(genes)}
    panel = {g: t for g, t in PANEL[cohort].items() if g in gidx and t in types}
    marg = np.array([types[int(P[:, gidx[g]].argmax())] for g in panel])
    mtgt = np.array(list(panel.values()))
    correct = marg == mtgt
    obs = int(correct.sum())
    pp, nmean, n95 = perm_p(marg, mtgt, obs, rng)
    per_lin = {t: f"{int(correct[mtgt == t].sum())}/{int((mtgt == t).sum())}" for t in types if (mtgt == t).any()}
    res = dict(n_genes_expressed=int(expressed.sum()), n_types=len(types),
               top_share=float(share[top]), top_type=types[top],
               share_by_type="; ".join(f"{t} {s:.3f}" for t, s in zip(types, share)),
               markers_correct=obs, markers_total=len(panel), perm_p=pp, null_mean=nmean, null_95=n95,
               per_lineage="; ".join(f"{k} {v}" for k, v in per_lin.items()))
    mtab = pd.DataFrame(dict(gene=list(panel), target=mtgt, argmax=marg, correct=correct.astype(int)))
    return res, mtab
