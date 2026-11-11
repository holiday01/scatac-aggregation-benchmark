"""
01: Marker concordance with marker genes taken from independent references (added 2026-11-11
in response to reviewer comment R1-2: the held-out panels were chosen by the author).

Markers come from CellTypist logistic-regression models trained on independent atlases
(already on the server; no new download):
  PBMC   Immune_All_Low      (Dominguez Conde et al., Science 2022)
  brain  Adult_Human_PrefrontalCortex (Ma et al., Science 2022)
  HCC    Healthy_Human_Liver (Guilliams et al., Cell 2022)
  COAD   Human_Colorectal_Cancer (Lee et al., Nat Genet 2020)
  NSCLC  Human_Lung_Atlas    (Sikkema et al., Nat Med 2023)
Each of our cell types is mapped to model classes by fixed name rules (written before any
result was seen and saved to marker_mapping.csv). For our type t and gene g the specificity
score is  mean coef over t's classes  minus  the largest mean coef over any other type's classes;
the top N genes with a positive score form t's panel (N = 10 primary, 20 sensitivity).
Excluded: genes used to assign the labels (labelling lists of PBMC and brain, the 72-gene COAD and
35-gene NSCLC scoring sets, every gene of the HCC label check), and genes absent from either
gene-activity matrix of the dataset, so both models are scored on the same panel.
Scoring is the paper's: a marker is correct when its per-gene maximum falls in its own type;
permutation null by shuffling the marker-to-type assignment over the observed maxima (10,000).
Out: independent_markers.csv (per dataset x model x strategy x N), marker_panels.csv, marker_mapping.csv
"""
import sys, re, time, os, warnings
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ADD = HERE.parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "tools"))
from benchmark_aggregation import m_raw, m_pseudobulk_cpm, m_depth_matched, m_cpm, m_lognorm, MATRICES, PANEL  # noqa
from celltypist import models as ctm
BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "/mnt/10t/scrna_atac"))
MODELS_DIR = Path(os.path.expanduser("~/.celltypist/data/models"))
SEED, NPERM = 42, 10_000
def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)
def P_of(r): return r[0] if isinstance(r, tuple) else r

# ── our type -> regex over model class names (case-insensitive, full match) ──
MAP = {
    "PBMC": ("Immune_All_Low.pkl", {
        "B": r"(B cells|Naive B cells|Memory B cells|Age-associated B cells)",
        "DC": r"(DC1|DC2|DC3)",
        "Mono_CD14": r"Classical monocytes",
        "Mono_CD16": r"Non-classical monocytes",
        "NK": r"(CD16\+ NK cells|CD16- NK cells|NK cells)",
        "T_CD4": r"(Tcm/Naive helper T cells|Tem/Effector helper T cells|Regulatory T cells|Follicular helper T cells)",
        "T_CD8": r"(Tcm/Naive cytotoxic T cells|Tem/Temra cytotoxic T cells|Tem/Trm cytotoxic T cells)",
        "pDC": r"pDC"}),
    "brain": ("Adult_Human_PrefrontalCortex.pkl", {
        "Astrocyte": r"Astro .*", "Excitatory_neuron": r"L\d.*", "Inhibitory_neuron": r"InN .*",
        "OPC": r"OPC .*", "Oligodendrocyte": r"Oligo .*"}),
    "HCC": ("Healthy_Human_Liver.pkl", {
        "B_cell": r"B cells", "CD4_T": r"T cells", "NK_cytotoxic_T": r"(Circulating NK/NKT|Resident NK)",
        "Endothelial_stromal": r"(Endothelial cells|Fibroblasts)", "Hepatocyte": r"Hepatocytes",
        "Macrophage": r"(Macrophages|Mono\+mono derived cells)"}),
    "COAD": ("Human_Colorectal_Cancer.pkl", {
        "B_cell": r"CD19\+CD20\+ B", "Endothelial": r"(Lymphatic ECs|Proliferative ECs|Stalk-like ECs|Tip-like ECs)",
        "Epithelial": r"(CMS1|CMS2|CMS3|CMS4|Goblet cells|Mature Enterocytes type 1|Mature Enterocytes type 2|Stem-like/TA|Intermediate)",
        "Fibroblast": r"(Myofibroblasts|Stromal 1|Stromal 2|Stromal 3)", "Macrophage": r"(SPP1\+|Pro-inflammatory)",
        "Mast": r"Mast cells", "Other_Myeloid": r"cDC",
        "T_NK": r"(CD4\+ T cells|CD8\+ T cells|NK cells|Regulatory T cells|T follicular helper cells|T helper 17 cells|gamma delta T cells)"}),
    "NSCLC": ("Human_Lung_Atlas.pkl", {
        "B_cell": r"B cells", "Endothelial": r"(EC .*|Lymphatic EC .*)",
        "Epithelial": r"(AT0|AT1|AT2|Basal resting|Suprabasal|Club \(non-nasal\)|Goblet \(bronchial\)|Goblet \(subsegmental\)|Multiciliated \(non-nasal\)|pre-TB secretory)",
        "Fibroblast": r"(Adventitial fibroblasts|Alveolar fibroblasts|Peribronchial fibroblasts|Subpleural fibroblasts|Myofibroblasts)",
        "Macrophage": r"(Alveolar macrophages|Alveolar Mph .*|Interstitial Mph perivascular|Monocyte-derived Mph)",
        "T_NK": r"(CD4 T cells|CD8 T cells|NK cells)"}),
}

# ── genes used to assign labels (excluded from the panels) ──
def lines(p): return set(Path(p).read_text().split())
COAD_LABEL = {"CD68", "MRC1", "CSF1R", "CD163", "MARCO", "C1QA", "C1QB", "C1QC", "APOE", "TREM2", "VSIG4", "STAB1", "MSR1",
              "CD14", "FCGR3A", "FCN1", "S100A8", "S100A9", "VCAN", "LYZ", "CSTA", "IL1B", "EREG", "CXCR2", "CXCR1", "FPR1",
              "FPR2", "CSF3R", "ELANE", "MPO", "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "NKG7", "KLRB1", "GZMA", "GZMB",
              "PRF1", "GNLY", "IL7R", "CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "JCHAIN", "IGHA1", "IGHG1", "TPSAB1", "TPSB2",
              "KIT", "CPA3", "EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "CDX2", "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
              "PDGFRA", "FAP", "PECAM1", "VWF", "CDH5", "KDR"}
NSCLC_LABEL = {"CD68", "CSF1R", "MRC1", "MSR1", "AIF1", "ITGAM", "MARCO", "MS4A1", "CD79A", "CD79B", "BANK1", "CD19", "PAX5",
               "CD3D", "CD3E", "CD2", "CD8A", "IL7R", "THEMIS", "NKG7", "GNLY", "KLRD1", "EPCAM", "KRT8", "KRT18", "CDH1",
               "NKX2-1", "SFTPC", "PECAM1", "VWF", "CDH5", "COL1A1", "COL1A2", "DCN", "LUM"}
assert len(COAD_LABEL) == 72 and len(NSCLC_LABEL) == 35
EXCLUDE = {"PBMC": lines(ADD / "multiome" / "labelling_genes.txt"),
           "brain": lines(ADD / "normal_tissues" / "brain3k_labelling_genes.txt"),
           "HCC": set(pd.read_csv(ROOT / "results" / "hcc_label_check.csv").gene),
           "COAD": COAD_LABEL, "NSCLC": NSCLC_LABEL}


def load_rna_labelled(lab_csv, paths):
    lab_df = pd.read_csv(lab_csv, index_col=0); lab_df = lab_df[lab_df.cell_type != "dropped"]
    out = []
    for model, p in paths:
        a = ad.read_h5ad(p); a = a[lab_df.index.intersection(a.obs_names)].copy()
        out.append((model, sp.csr_matrix(a.X).astype(np.float32), a.var_names.astype(str).values,
                    lab_df.loc[a.obs_names, "cell_type"].values.astype(str)))
    return out


def load_tumor(cohort):
    out = []
    for c, model, rel, labcol in MATRICES:
        if c != cohort or model not in ("proximal", "enhancer"): continue
        a = ad.read_h5ad(BASE / rel)
        lab_all = a.obs[labcol].astype(str).values
        keep = ~pd.isna(a.obs[labcol].values) & (lab_all != "nan")
        a = a[keep].copy(); lab = a.obs[labcol].astype(str).values
        if cohort == "HCC":   # same correction as the benchmark (tools/benchmark_aggregation.py)
            lab = np.where(lab == "B_cell", "Endothelial_stromal", np.where(lab == "DC", "B_cell", lab))
        out.append(("count" if model == "proximal" else "weighted", sp.csr_matrix(a.X).astype(np.float32),
                    a.var_names.astype(str).values, lab))
    return out


LOADERS = {
    "PBMC": lambda: load_rna_labelled(ADD / "multiome" / "rna_labels.csv",
        [("count", BASE / "data/processed/multiome_pbmc10k/pbmc10k_geneactivity_proximal.h5ad"),
         ("weighted", BASE / "data/processed/multiome_pbmc10k/pbmc10k_geneactivity_enhancer.h5ad")]),
    "brain": lambda: load_rna_labelled(ADD / "normal_tissues" / "brain3k_rna_labels.csv",
        [("count", BASE / "data/processed/multiome_brain3k/brain3k_geneactivity_proximal.h5ad"),
         ("weighted", BASE / "data/processed/multiome_brain3k/brain3k_geneactivity_enhancer.h5ad")]),
    "HCC": lambda: load_tumor("HCC"), "COAD": lambda: load_tumor("COAD"), "NSCLC": lambda: load_tumor("LUAD"),
}

rng = np.random.default_rng(SEED)
res, panels_out, map_out = [], [], []
for ds, (mfile, rules) in MAP.items():
    mats = LOADERS[ds]()
    types = sorted(set(mats[0][3]))
    assert set(rules) == set(types), (ds, sorted(rules), types)
    m = ctm.Model.load(str(MODELS_DIR / mfile))
    coef = m.classifier.coef_; feats = np.asarray(m.features).astype(str); classes = list(m.cell_types)
    groups = {}
    for t, rx in rules.items():
        idx = [i for i, c in enumerate(classes) if re.fullmatch(rx, c, flags=re.I)]
        assert idx, (ds, t, rx)
        groups[t] = idx
        map_out.append(dict(dataset=ds, model_file=mfile, our_type=t, model_classes="; ".join(classes[i] for i in idx)))
    used = set(i for v in groups.values() for i in v)
    assert len(used) == sum(len(v) for v in groups.values()), (ds, "a model class maps to two types")
    G = np.vstack([coef[groups[t]].mean(0) for t in types])            # types x features
    common = set(mats[0][2])
    for _, X, g, _ in mats[1:]: common &= set(g)
    expressed = set.intersection(*[set(g[np.asarray(X.sum(0)).ravel() > 0]) for _, X, g, _ in mats])
    eligible = np.array([(f in common) and (f in expressed) and (f not in EXCLUDE[ds]) for f in feats])
    author = set(PANEL.get({"NSCLC": "LUAD"}.get(ds, ds), {})) if isinstance(PANEL, dict) else set()
    for N in (10, 20):
        panel = {}
        for k, t in enumerate(types):
            other = np.delete(G, k, axis=0).max(0)
            score = np.where(eligible, G[k] - other, -np.inf)
            top = [i for i in np.argsort(-score)[:N] if score[i] > 0]
            for i in top:
                panel[feats[i]] = t
                panels_out.append(dict(dataset=ds, N=N, our_type=t, gene=feats[i], score=float(score[i]),
                                       in_author_panel=feats[i] in author))
        for model, X, genes, lab in mats:
            gidx = {x: i for i, x in enumerate(genes)}
            mk = [x for x in panel if x in gidx]
            strategies = {"raw_mean": m_raw, "lognorm_mean": m_lognorm, "pseudobulk_cpm": m_pseudobulk_cpm, "cpm_mean": m_cpm}
            Ps = {s: P_of(f(X, lab, types)) for s, f in strategies.items()}
            if model == "count":
                Ps["depth_matched"] = P_of(m_depth_matched(X, lab, types, rng))
            for s, P in Ps.items():
                arg = np.array([types[int(P[:, gidx[x]].argmax())] for x in mk])
                tgt = np.array([panel[x] for x in mk])
                obs = int((arg == tgt).sum())
                null = np.array([(rng.permutation(arg) == tgt).sum() for _ in range(NPERM)])
                per = "; ".join(f"{t} {int(((arg == tgt) & (tgt == t)).sum())}/{int((tgt == t).sum())}" for t in types)
                res.append(dict(dataset=ds, model=model, N=N, strategy=s, correct=obs, total=len(mk),
                                null_mean=float(null.mean()), perm_p=float((null >= obs).mean()), per_type=per,
                                n_in_author_panel=int(sum(x in author for x in mk))))
            log(f"{ds} {model} N={N}: " + ", ".join(f"{r['strategy']} {r['correct']}/{r['total']}"
                                                    for r in res if r['dataset'] == ds and r['model'] == model and r['N'] == N))
    pd.DataFrame(res).to_csv(HERE / "independent_markers.csv", index=False)
    pd.DataFrame(panels_out).to_csv(HERE / "marker_panels.csv", index=False)
    pd.DataFrame(map_out).to_csv(HERE / "marker_mapping.csv", index=False)
log("done")
