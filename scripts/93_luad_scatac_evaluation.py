"""
Script 93: Is LUAD (GSE270148) a valid replacement for the broken COAD arm?
==========================================================================
COAD must be dropped: 0/1434 macrophages and 0/1759 T_NK cells survive a
>=1000-fragment QC (Script 92). We need a third cancer with REAL cells,
all three lineages (Macrophage / B / CD8-T), a TCGA survival cohort and a
scRNA reference.

LUAD GSE270148 is the leading candidate (peak matrices + TCGA_LUAD + scRNA_LUAD
already local). But the manuscript rejected it on the claim of a
"B-cell library-prep artefact (4 of 9 patients with B-cell fractions >70%)".
That claim was produced by the SAME annotation machinery that we now know is
broken, so it is re-tested here from scratch:

  * gene activity computed directly from each sample's peak matrix
    (peaks overlapping gene body + 2 kb upstream — the Signac model),
    using the peak coordinates stored in the 10x h5 (no fragments needed)
  * cells annotated by argmax of canonical lineage marker scores
  * per-sample composition + depth QC reported

NOTE two data defects found in GSE270148 itself:
  - GSM8069382 (Lu931) and GSM8069383 (Lu934) are BYTE-IDENTICAL on GEO
    (same md5) -> a duplicate submission. Only 8 unique samples exist.
  - GSM8069377 (Lu883) has median 304 counts/cell (7% pass >=1000) -> unusable.
  => 7 usable samples.

Output: results/validation/luad_scatac_evaluation.csv
"""
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
from pathlib import Path
from datetime import datetime
import glob, os, warnings
warnings.filterwarnings('ignore')

BASE = Path(os.environ.get("SCRNA_ATAC_BASE", "."))
H5   = sorted(glob.glob(str(BASE / "data/raw/scATAC_LUAD/GSE270148/*.h5")))
GTF  = os.environ.get("GENCODE_GTF", "data/reference/gencode.v38.annotation.gtf")
OUT  = BASE / "results/validation"
OUT.mkdir(parents=True, exist_ok=True)

MIN_COUNTS = 1000          # standard scATAC cell QC
DUPLICATE  = 'GSM8069383'  # byte-identical to GSM8069382 on GEO
LOWQC      = 'GSM8069377'  # median 304 counts/cell

MARKERS = {
    'Macrophage': ['CD68', 'CSF1R', 'MRC1', 'MSR1', 'AIF1', 'ITGAM', 'MARCO', 'MSR1'],
    'B_cell':     ['MS4A1', 'CD79A', 'CD79B', 'BANK1', 'CD19', 'PAX5'],
    'T_cell':     ['CD3D', 'CD3E', 'CD2', 'CD8A', 'IL7R', 'THEMIS'],
    'NK':         ['NKG7', 'GNLY', 'KLRD1', 'PRF1'],
    'Epithelial': ['EPCAM', 'KRT8', 'KRT18', 'CDH1', 'NKX2-1', 'SFTPC'],
    'Endothelial':['PECAM1', 'VWF', 'CDH5'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
}


def log(m):
    print(f"[{datetime.now():%H:%M:%S}] {m}", flush=True)


# ── gene bodies (+2 kb upstream), protein-coding ─────────────────────────────
log("Parsing GENCODE v38 gene bodies")
rows = []
with open(GTF) as f:
    for line in f:
        if line[0] == '#':
            continue
        p = line.split('\t')
        if p[2] != 'gene' or 'gene_type "protein_coding"' not in p[8]:
            continue
        chrom, s, e, strand = p[0], int(p[3]), int(p[4]), p[6]
        s2, e2 = (s - 2000, e) if strand == '+' else (s, e + 2000)   # +2kb upstream
        rows.append((chrom, s2, e2, p[8].split('gene_name "')[1].split('"')[0]))
genes = pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'gene']).drop_duplicates('gene')
want = sorted({g for v in MARKERS.values() for g in v})
genes = genes[genes.gene.isin(want)].reset_index(drop=True)
log(f"  marker genes located: {len(genes)}/{len(want)}")

records = []
for f in H5:
    gsm = os.path.basename(f).split('_')[0]
    tag = ''
    if gsm == DUPLICATE:
        log(f"{gsm}: SKIP — byte-identical duplicate of GSM8069382 on GEO"); continue
    a = sc.read_10x_h5(f, gex_only=False)
    a.var_names_make_unique()

    depth = np.asarray(a.X.sum(1)).ravel()
    keep = depth >= MIN_COUNTS
    if keep.sum() < 100:
        log(f"{gsm}: SKIP — only {keep.sum()} cells pass >={MIN_COUNTS} counts"); continue
    a = a[keep].copy()

    # peak coords
    v = a.var.copy()
    if {'chrom', 'chromStart', 'chromEnd'} <= set(v.columns):
        pc = v[['chrom', 'chromStart', 'chromEnd']].copy()
        pc.columns = ['chrom', 'start', 'end']
    else:                                    # 'chr1:100-200' style var_names
        n = pd.Series(a.var_names)
        pc = pd.DataFrame({'chrom': n.str.split(':').str[0]})
        se = n.str.split(':').str[1].str.split('-')
        pc['start'] = se.str[0].astype(int); pc['end'] = se.str[1].astype(int)
    pc = pc.reset_index(drop=True)
    pc['start'] = pc.start.astype(int); pc['end'] = pc.end.astype(int)
    pc['idx'] = np.arange(len(pc))

    # peak x marker-gene overlap weight matrix
    wi, wj = [], []
    by = {c: g for c, g in pc.groupby('chrom')}
    for gi, g in enumerate(genes.itertuples(index=False)):
        ch = by.get(g.chrom)
        if ch is None:
            continue
        m = ch[(ch.end >= g.start) & (ch.start <= g.end)]     # overlap
        wi.extend(m.idx.values); wj.extend([gi] * len(m))
    W = sp.csr_matrix((np.ones(len(wi)), (wi, wj)), shape=(len(pc), len(genes)))
    G = np.asarray((a.X @ W).todense())                        # cells x marker genes
    G = pd.DataFrame(G, columns=genes.gene.values)

    # depth-normalise, then z-score each gene, then lineage score = mean z
    d = np.asarray(a.X.sum(1)).ravel()
    Gn = G.div(d, axis=0) * 1e4
    Gz = (Gn - Gn.mean()) / (Gn.std() + 1e-9)
    scores = pd.DataFrame({
        lin: Gz[[g for g in gl if g in Gz.columns]].mean(axis=1)
        for lin, gl in MARKERS.items()
    })
    call = scores.idxmax(axis=1)
    # merge NK into T_cell for the lineage tally (paper uses CD8/NK jointly)
    call = call.replace({'NK': 'T_cell'})

    comp = call.value_counts(normalize=True)
    n = len(call)
    records.append(dict(
        gsm=gsm, sample=os.path.basename(f).split('_')[4],
        cells_pass_qc=n, median_counts=float(np.median(d)),
        pct_Macrophage=round(100 * comp.get('Macrophage', 0), 1),
        pct_B_cell=round(100 * comp.get('B_cell', 0), 1),
        pct_T_NK=round(100 * comp.get('T_cell', 0), 1),
        pct_Epithelial=round(100 * comp.get('Epithelial', 0), 1),
        n_Macrophage=int((call == 'Macrophage').sum()),
        n_B_cell=int((call == 'B_cell').sum()),
        n_T_NK=int((call == 'T_cell').sum()),
    ))
    log(f"  {gsm} {records[-1]['sample']:6} n={n:6d} med={np.median(d):7.0f}  "
        f"Mac={records[-1]['pct_Macrophage']:5.1f}%  B={records[-1]['pct_B_cell']:5.1f}%  "
        f"T/NK={records[-1]['pct_T_NK']:5.1f}%  Epi={records[-1]['pct_Epithelial']:5.1f}%")

df = pd.DataFrame(records)
df.to_csv(OUT / "luad_scatac_evaluation.csv", index=False)

print("\n" + "=" * 78)
print("LUAD GSE270148 — usable samples after removing GEO duplicate + low-QC sample")
print("=" * 78)
print(df.to_string(index=False))
print()
print(f"  usable samples          : {len(df)}")
print(f"  total cells (>={MIN_COUNTS} counts): {df.cells_pass_qc.sum():,}")
print(f"  TOTAL Macrophage cells  : {df.n_Macrophage.sum():,}")
print(f"  TOTAL B cells           : {df.n_B_cell.sum():,}")
print(f"  TOTAL T/NK cells        : {df.n_T_NK.sum():,}")
print()
print(f"  patients with B-cell fraction >70%: {(df.pct_B_cell > 70).sum()} / {len(df)}"
      f"   <- manuscript claimed 4 of 9")
log("Script 93 done.")
