"""
01: Gene-activity matrices for the adult non-tumor tissues of the human chromatin
    accessibility atlas (Zhang et al. 2021, Cell; GSE184462), built from the deposited
    per-sample fragment files with the same two models used for the tumor datasets and
    the PBMC multiome (additions/multiome/01_build_gene_activity.py):

  proximal  : snapatac2 make_gene_matrix(upstream=2000, downstream=0, include_gene_body=True)
  enhancer  : distance-weighted gene score over 5 kb tiles,
              w = exp(-d/5kb) within 100 kb of the gene body, full weight inside the body

Cells and labels are the atlas's own (GSE184462_metadata.tsv.gz, "cell type"); only
barcodes present in that table are imported. True per-cell depth = n_fragment recorded
by snapatac2 at import (never a matrix row sum).

Run:  python3 01_build_catlas.py liver lung            # one or more organs
      python3 01_build_catlas.py --list                # organs with all fragment files on disk
Out:  data/processed/catlas_adult/{organ}_geneactivity_{proximal,enhancer}.h5ad
"""
import sys, time, re, warnings
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, hdf5plugin  # noqa
import snapatac2 as snap
warnings.filterwarnings("ignore")

RAW  = Path("/mnt/10t/scrna_atac/data/raw/catlas_adult")
FRAG = RAW / "fragments"
GTF  = Path("/mnt/10t/holiday/hnsc_analysis/gencode.v38.annotation.gtf")
OUT  = Path("/mnt/10t/scrna_atac/data/processed/catlas_adult"); OUT.mkdir(parents=True, exist_ok=True)
TMP  = OUT / "_per_sample"; TMP.mkdir(exist_ok=True)
BIN, EXTEND, DECAY = 5000, 100_000, 5000
N_JOBS = 6


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def load_meta():
    m = pd.read_csv(RAW / "GSE184462_metadata.tsv.gz", sep="\t")
    m = m[m["Life stage"] == "Adult"].copy()
    m["organ"] = m["tissue"].str.replace(r"_SM-\w+$", "", regex=True)
    # the atlas calls the two frontal-cortex replicates of donor UMB4540 Human_brain_1 / Human_brain_2
    # (checked by barcode overlap with GSM*_UMB4540_snATAC_frontal_cortex_rep1 / rep2)
    m.loc[m["tissue"].isin(BRAIN), "organ"] = "frontal_cortex"
    m["barcode"] = m["cellID"].str.split("+").str[1]
    m["donor"] = m["tissue"].str.extract(r"_(SM-\w+)$")[0]
    m.loc[m["tissue"].isin(BRAIN), "donor"] = "UMB4540"
    return m


BRAIN = {"Human_brain_1": "UMB4540_snATAC_frontal_cortex_rep1", "Human_brain_2": "UMB4540_snATAC_frontal_cortex_rep2"}


def frag_file(sample):
    """metadata sample 'liver_SM-A8WNZ_1' -> GSM*_liver_SM-A8WNZ_rep1_fragments.bed.gz"""
    if sample.rsplit("_", 1)[0] in BRAIN:
        hits = list(FRAG.glob(f"GSM*_{BRAIN[sample.rsplit('_', 1)[0]]}_fragments.bed.gz"))
        return hits[0] if len(hits) == 1 else None
    stem, rep = sample.rsplit("_", 1)
    hits = list(FRAG.glob(f"GSM*_{stem}_rep{rep}_fragments.bed.gz"))
    return hits[0] if len(hits) == 1 else None


def gene_table(chroms):
    rows = []
    with open(GTF) as f:
        for line in f:
            if line[0] == "#": continue
            p = line.split("\t")
            if p[2] != "gene" or 'gene_type "protein_coding"' not in p[8]: continue
            rows.append((p[0], int(p[3]), int(p[4]), p[6], p[8].split('gene_name "')[1].split('"')[0]))
    genes = pd.DataFrame(rows, columns=["chrom", "start", "end", "strand", "gene"]).drop_duplicates("gene")
    return genes[genes.chrom.isin(chroms)].reset_index(drop=True)


def tile_weights(tile_names):
    """tiles x genes weight matrix, identical to the PBMC multiome build; cached."""
    cache = OUT / "_W_tiles5kb.npz"
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        if list(z["tiles"]) == list(tile_names):
            return sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])), list(z["genes"])
    tc = pd.DataFrame({"name": tile_names})
    tc["chrom"] = tc.name.str.split(":").str[0]
    se = tc.name.str.split(":").str[1].str.split("-")
    tc["start"] = se.str[0].astype(int); tc["end"] = se.str[1].astype(int)
    tc["mid"] = (tc.start + tc.end) // 2; tc["idx"] = np.arange(len(tc))
    genes = gene_table(tc.chrom.unique())
    by_chrom = {c: g.sort_values("mid") for c, g in tc.groupby("chrom")}
    wi, wj, wv = [], [], []
    for gi, g in enumerate(genes.itertuples(index=False)):
        ch = by_chrom.get(g.chrom)
        if ch is None: continue
        m = ch[(ch.mid >= g.start - EXTEND) & (ch.mid <= g.end + EXTEND)]
        if len(m) == 0: continue
        d = np.where(m.mid.values < g.start, g.start - m.mid.values,
                     np.where(m.mid.values > g.end, m.mid.values - g.end, 0))
        wi.extend(m.idx.values); wj.extend([gi] * len(m)); wv.extend(np.exp(-d / DECAY))
    W = sp.csr_matrix((wv, (wi, wj)), shape=(len(tc), len(genes)))
    np.savez(cache, data=W.data, indices=W.indices, indptr=W.indptr, shape=np.array(W.shape),
             tiles=np.array(tile_names, dtype=object), genes=genes.gene.values.astype(object))
    return W, list(genes.gene.values)


def build_sample(sample, meta_s):
    """-> (proximal AnnData, enhancer AnnData) for one sample, cached on disk."""
    pp, pe = TMP / f"{sample}_proximal.h5ad", TMP / f"{sample}_enhancer.h5ad"
    if pp.exists() and pe.exists():
        return ad.read_h5ad(pp), ad.read_h5ad(pe)
    f = frag_file(sample)
    barcodes = list(meta_s.barcode)
    fh5 = TMP / f"{sample}_fragments.h5ad"
    if fh5.exists(): fh5.unlink()
    data = snap.pp.import_fragments(f, chrom_sizes=snap.genome.hg38, file=str(fh5), whitelist=barcodes,
                                    min_num_fragments=0, sorted_by_barcode=False, n_jobs=N_JOBS)
    obs = list(data.obs_names)
    nfrag = pd.Series(np.asarray(data.obs["n_fragment"][:]), index=obs)
    log(f"  {sample}: imported {len(obs)}/{len(barcodes)} labelled barcodes; n_fragment median {nfrag.median():.0f}")

    gm = snap.pp.make_gene_matrix(data, gene_anno=GTF, upstream=2000, downstream=0,
                                  include_gene_body=True, file=str(TMP / f"{sample}_gm.h5ad"))
    P = ad.AnnData(X=sp.csr_matrix(gm.X[:]).astype(np.float32),
                   obs=pd.DataFrame(index=list(gm.obs_names)), var=pd.DataFrame(index=list(gm.var_names)))
    gm.close(); (TMP / f"{sample}_gm.h5ad").unlink()
    P.var_names_make_unique()

    snap.pp.add_tile_matrix(data, bin_size=BIN, exclude_chroms=["chrM", "chrY", "M", "Y"], n_jobs=N_JOBS)
    T = sp.csr_matrix(data.X[:]); tiles = list(data.var_names)
    data.close(); fh5.unlink()
    W, gnames = tile_weights(tiles)
    E = ad.AnnData(X=(T @ W).tocsr().astype(np.float32), obs=pd.DataFrame(index=obs),
                   var=pd.DataFrame(index=gnames))
    E.var_names_make_unique()

    ms = meta_s.set_index("barcode")
    for a in (P, E):
        a.obs["n_fragment"] = nfrag.reindex(a.obs_names).values
        a.obs["cell_type"] = ms["cell type"].reindex(a.obs_names).values
        a.obs["sample"] = sample
        a.obs["donor"] = ms["donor"].reindex(a.obs_names).values
        a.obs["atlas_logUMI"] = ms["logUMI"].reindex(a.obs_names).values
        a.obs_names = [f"{sample}+{b}" for b in a.obs_names]
    P.write_h5ad(pp); E.write_h5ad(pe)
    return P, E


def main():
    meta = load_meta()
    if "--list" in sys.argv:
        for organ, m in meta.groupby("organ"):
            samples = sorted(m["sample"].unique())
            have = [s for s in samples if frag_file(s) is not None]
            if have: print(f"{organ}\t{len(have)}/{len(samples)} samples on disk\t{len(m)} cells")
        return
    for organ in sys.argv[1:]:
        m = meta[meta.organ == organ]
        samples = [s for s in sorted(m["sample"].unique()) if frag_file(s) is not None]
        if not samples:
            log(f"{organ}: no fragment files on disk, skipped"); continue
        log(f"{organ}: {len(samples)} samples, {len(m)} labelled cells")
        parts = [build_sample(s, m[m["sample"] == s]) for s in samples]
        for k, model in enumerate(("proximal", "enhancer")):
            a = ad.concat([p[k] for p in parts], join="inner")
            a.write_h5ad(OUT / f"{organ}_geneactivity_{model}.h5ad")
            log(f"  saved {organ} {model} {a.shape}")
    log("done")


if __name__ == "__main__":
    main()
