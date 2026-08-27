"""
Script 90: Rebuild HCC scATAC gene activity from fragments — correctly named,
           with a proximal (Signac-equivalent) and an enhancer-aware (distal) model.
=============================================================================
MOTIVATION (two separate problems this script resolves)

(1) ANNOTATION BUG. The canonical prior `atac_hcc_geneactivity.h5ad` was built by
    `scripts/preprocessing/04_geneactivity_names.py`, which took the GEO-supplied
    GSE227265_GeneActivity.mtx (no row names) and *guessed* gene names by sorting
    GENCODE v27 protein-coding genes by chr+start ("ArchR convention").
    But Craig et al. (Cell Rep 2023, GSE227265) actually produced that matrix with
    Signac::GeneActivity() using EnsDb.Hsapiens.v86 — a different annotation and a
    different row order. The guess is wrong: canonical lineage markers land in the
    wrong cell type (ALB argmax = DC; CSF1R argmax = Hepatocyte; only 5/22 correct
    once cells are depth-normalised before the cell-type mean -- see
    apbc2026/tools/audit_marker_concordance.py. Un-normalised the count reads 6/22;
    that variant is superseded, and the difference is documented in
    apbc2026/STATUS_NUMBERS_SUPERSEDED.md).

(2) PROXIMAL-ONLY PRIOR. Signac GeneActivity counts gene body + 2 kb upstream, so it
    cannot see distal enhancers. The manuscript's mechanism claim (the prior demotes
    enhancer-driven markers such as CD163/CD79B) rests on this, and its Limitations
    concede the null may be an artefact of the proxy. Testing that requires an
    enhancer-aware prior, which has never been built.
    NOTE (2026-08-27): that mechanism claim did not survive. Depth-normalised,
    CD163 is correctly assigned under BOTH priors and the whole macrophage panel
    goes 1/6 -> 6/6; only CD79B and MS4A1 remain wrong. The demotion story was an
    artefact of averaging un-normalised cells.

This script rebuilds gene activity DIRECTLY FROM FRAGMENTS with real coordinates:

  A. PROXIMAL  — snapatac2 make_gene_matrix(upstream=2000, downstream=0,
                 include_gene_body=True) == Signac GeneActivity model, correctly named.
                 Isolates the naming bug: same model as GEO, right names.
  B. ENHANCER  — ArchR-style distance-weighted gene score over 5 kb tiles:
                 w(tile,gene) = exp(-d / DECAY) for tiles within EXTEND of the gene
                 body, plus full weight inside the gene body. Sees distal peaks.

Comparing A vs B isolates the *enhancer* question with names held correct.
Comparing OLD vs A isolates the *naming* question with the model held constant.

Inputs : GSE227265_fragments_AllSamples.tsv.gz (coordinate-sorted, bgzipped)
         atac_hcc_geneactivity.h5ad  (for the 12,029 HCC cell barcodes + cell_type)
         gencode.v38.annotation.gtf  (hg38)
Outputs: data/processed/integrated/atac_hcc_geneactivity_proximal.h5ad
         data/processed/integrated/atac_hcc_geneactivity_enhancer.h5ad
         results/validation/gene_activity_rebuild_validation.csv
"""
import numpy as np, pandas as pd, anndata as ad
import scipy.sparse as sp
import snapatac2 as snap
from pathlib import Path
from datetime import datetime
import sys, warnings
warnings.filterwarnings('ignore')

BASE     = Path("/mnt/10t/scrna_atac")
FRAG     = Path("/mnt/10t/assi_result/GEO_DATA/GSE227265/suppl/GSE227265_fragments_AllSamples.tsv.gz")
GTF      = Path("/mnt/10t/holiday/hnsc_analysis/gencode.v38.annotation.gtf")
OLD_H5   = BASE / "data/processed/integrated/atac_hcc_geneactivity.h5ad"
OUT_DIR  = BASE / "data/processed/integrated"
VAL_DIR  = BASE / "results/validation"
WORK     = OUT_DIR / "_s90_work"
WORK.mkdir(parents=True, exist_ok=True)
VAL_DIR.mkdir(parents=True, exist_ok=True)

FRAG_H5  = WORK / "hcc_fragments.h5ad"
TILE_H5  = WORK / "hcc_tiles.h5ad"

# ArchR-style gene-score parameters
BIN       = 5000      # tile size (bp)
EXTEND    = 100_000   # consider tiles within 100 kb of the gene body
DECAY     = 5000      # exponential decay constant (bp)

STAGE = sys.argv[1] if len(sys.argv) > 1 else 'all'


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


# ── cells ────────────────────────────────────────────────────────────────────
log("Loading HCC cell whitelist from canonical h5ad")
old = ad.read_h5ad(OLD_H5)
barcodes = list(old.obs_names)          # e.g. 'TTGTTCAAGGAATGGA-12'
cell_type = old.obs['cell_type'].astype(str)
log(f"  {len(barcodes)} HCC cells, {cell_type.nunique()} cell types")


# ── A. import fragments ──────────────────────────────────────────────────────
if STAGE in ('all', 'import') and not FRAG_H5.exists():
    log(f"Importing fragments (coordinate-sorted, NOT barcode-sorted) -> {FRAG_H5}")
    data = snap.pp.import_fragments(
        FRAG,
        chrom_sizes=snap.genome.hg38,
        file=str(FRAG_H5),
        whitelist=barcodes,
        min_num_fragments=0,        # keep every whitelisted cell
        sorted_by_barcode=False,    # this file is sorted by coordinate
        n_jobs=8,
    )
    log(f"  imported: {data.shape}")
    data.close()
else:
    log(f"  fragments h5ad exists: {FRAG_H5.exists()}")


# ── B. proximal gene matrix (Signac-equivalent, correctly named) ─────────────
if STAGE in ('all', 'proximal'):
    log("Building PROXIMAL gene matrix (gene body + 2 kb upstream) == Signac model")
    data = snap.read(str(FRAG_H5), backed='r')
    gm = snap.pp.make_gene_matrix(
        data, gene_anno=GTF,
        upstream=2000, downstream=0, include_gene_body=True,
        file=str(WORK / "gm_proximal.h5ad"),
    )
    X = gm.X[:] if not isinstance(gm.X, np.ndarray) else gm.X
    obs_names = list(gm.obs_names)
    var_names = list(gm.var_names)
    gm.close(); data.close()

    a = ad.AnnData(X=sp.csr_matrix(X),
                   obs=pd.DataFrame(index=obs_names),
                   var=pd.DataFrame(index=var_names))
    a.var['gene_name'] = a.var_names
    a = a[[b for b in barcodes if b in a.obs_names]].copy()
    a.obs['cell_type'] = cell_type.reindex(a.obs_names).values
    a.obs['sample_id'] = [b.rsplit('-', 1)[-1] for b in a.obs_names]
    a.obs['barcode']   = [b.rsplit('-', 1)[0]  for b in a.obs_names]
    a.var_names_make_unique()
    a.write_h5ad(OUT_DIR / "atac_hcc_geneactivity_proximal.h5ad")
    log(f"  saved proximal: {a.shape}")


# ── C. enhancer-aware gene matrix (ArchR-style distance weighting) ───────────
if STAGE in ('all', 'enhancer'):
    log(f"Building tile matrix (bin={BIN}) for enhancer-aware score")
    # add_tile_matrix(inplace=True) writes X into the fragments h5ad itself
    # (the `file=` arg only applies when inplace=False), so open r+ and read back.
    data = snap.read(str(FRAG_H5), backed='r+')
    if data.X is None or data.shape[1] == 0:
        snap.pp.add_tile_matrix(data, bin_size=BIN,
                                exclude_chroms=['chrM', 'chrY', 'M', 'Y'], n_jobs=8)
    T = sp.csr_matrix(data.X[:])        # cells x tiles
    tile_names = list(data.var_names)   # 'chr1:0-5000'
    tobs = list(data.obs_names)
    data.close()
    log(f"  tile matrix: {T.shape}")

    # tile coordinates
    tc = pd.DataFrame({'name': tile_names})
    tc['chrom'] = tc.name.str.split(':').str[0]
    se = tc.name.str.split(':').str[1].str.split('-')
    tc['start'] = se.str[0].astype(int)
    tc['end']   = se.str[1].astype(int)
    tc['mid']   = (tc.start + tc.end) // 2
    tc['idx']   = np.arange(len(tc))

    # gene bodies from GTF (protein-coding)
    log("Parsing gene bodies from GENCODE v38")
    rows = []
    with open(GTF) as f:
        for line in f:
            if line[0] == '#':
                continue
            p = line.split('\t')
            if p[2] != 'gene':
                continue
            attr = p[8]
            if 'gene_type "protein_coding"' not in attr:
                continue
            gname = attr.split('gene_name "')[1].split('"')[0]
            rows.append((p[0], int(p[3]), int(p[4]), p[6], gname))
    genes = pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'strand', 'gene'])
    genes = genes.drop_duplicates('gene')
    genes = genes[genes.chrom.isin(tc.chrom.unique())]
    log(f"  protein-coding genes: {len(genes)}")

    # weight matrix W (tiles x genes): exp(-d/DECAY), d=0 inside gene body
    log(f"Building distance-weighted W (extend={EXTEND}, decay={DECAY})")
    by_chrom = {c: g.sort_values('mid') for c, g in tc.groupby('chrom')}
    wi, wj, wv = [], [], []
    for gi, g in enumerate(genes.itertuples(index=False)):
        ch = by_chrom.get(g.chrom)
        if ch is None:
            continue
        lo, hi = g.start - EXTEND, g.end + EXTEND
        m = ch[(ch.mid >= lo) & (ch.mid <= hi)]
        if len(m) == 0:
            continue
        d = np.where(m.mid.values < g.start, g.start - m.mid.values,
            np.where(m.mid.values > g.end,   m.mid.values - g.end, 0))
        w = np.exp(-d / DECAY)
        wi.extend(m.idx.values); wj.extend([gi] * len(m)); wv.extend(w)
    W = sp.csr_matrix((wv, (wi, wj)), shape=(len(tc), len(genes)))
    log(f"  W: {W.shape}, nnz={W.nnz:,}")

    log("GeneScore = Tiles @ W")
    G = (T @ W).tocsr()

    a = ad.AnnData(X=G,
                   obs=pd.DataFrame(index=tobs),
                   var=pd.DataFrame(index=genes.gene.values))
    a.var['gene_name'] = a.var_names
    a = a[[b for b in barcodes if b in a.obs_names]].copy()
    a.obs['cell_type'] = cell_type.reindex(a.obs_names).values
    a.obs['sample_id'] = [b.rsplit('-', 1)[-1] for b in a.obs_names]
    a.obs['barcode']   = [b.rsplit('-', 1)[0]  for b in a.obs_names]
    a.var_names_make_unique()
    a.write_h5ad(OUT_DIR / "atac_hcc_geneactivity_enhancer.h5ad")
    log(f"  saved enhancer: {a.shape}")


# ── D. validation: marker concordance across the three priors ────────────────
if STAGE in ('all', 'validate'):
    log("Validating: does each prior put canonical markers in the right lineage?")
    LIN = {'CD3D':'T','CD3E':'T','CD8A':'T','IL7R':'T','NKG7':'T','GNLY':'T',
           'MS4A1':'B','CD79A':'B','CD79B':'B','CD19':'B','BANK1':'B',
           'CSF1R':'Mac','CD68':'Mac','MRC1':'Mac','MSR1':'Mac','CD163':'Mac','ITGAM':'Mac',
           'ALB':'Hep','APOB':'Hep','HNF4A':'Hep','TTR':'Hep','APOA1':'Hep'}
    EXP = {'T':'NK_cytotoxic_T','B':'B_cell','Mac':'Macrophage','Hep':'Hepatocyte'}

    def concordance(h5, tag):
        a = ad.read_h5ad(h5)
        X = a.X.tocsc() if sp.issparse(a.X) else a.X
        ct = a.obs['cell_type'].astype(str).values
        out = []
        for g, lin in LIN.items():
            if g not in a.var_names:
                out.append((tag, g, lin, 'ABSENT', 0)); continue
            v = X[:, a.var_names.get_loc(g)]
            v = np.asarray(v.todense()).ravel() if sp.issparse(v) else np.asarray(v).ravel()
            m = pd.Series(v).groupby(ct).mean()
            am = m.idxmax()
            out.append((tag, g, lin, am, int(am == EXP[lin])))
        return out

    recs = []
    recs += concordance(OLD_H5, 'old_guessed_names')
    recs += concordance(OUT_DIR / "atac_hcc_geneactivity_proximal.h5ad", 'new_proximal')
    recs += concordance(OUT_DIR / "atac_hcc_geneactivity_enhancer.h5ad", 'new_enhancer')
    df = pd.DataFrame(recs, columns=['prior', 'gene', 'lineage', 'atac_argmax', 'correct'])
    df.to_csv(VAL_DIR / "gene_activity_rebuild_validation.csv", index=False)

    print()
    print(df.groupby('prior').correct.agg(['sum', 'count']).to_string())
    print()
    log(f"saved {VAL_DIR / 'gene_activity_rebuild_validation.csv'}")

log("Script 90 done.")
