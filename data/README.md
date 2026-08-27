# Data

No data files are stored in this repository. Everything analysed is public, and
the derived tables the manuscript reports are in [`../results/`](../results/).
Checksums of the deposited input files are in [`checksums.md`](checksums.md).

The gene-activity matrices are 0.6–1.6 GB each and exceed GitHub's file limit;
they are regenerated from the deposited files by
`../scripts/90_enhancer_aware_gene_activity.py` (HCC), `../scripts/92_coad_gene_activity_rebuild.py`
(colorectal) and `../scripts/94_luad_prior_and_reference.py` (lung).

Scripts resolve every data path from `$SCRNA_ATAC_BASE` (default: the repository root).
The GENCODE GTF is resolved from `$GENCODE_GTF` and the GSE227265 fragment file
from `$GSE227265_FRAGMENTS`.
Set it to the root you populate below.

## Accessions

| Cohort | Accession | Used for |
|---|---|---|
| HCC scATAC | [GSE227265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227265) | fragments (count and weighted matrices, true per-cell coverage); the deposited `GeneActivity.mtx` is the base for the constructed negative control (see below) |
| HCC scRNA | [GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614) | transcriptional reference |
| Colorectal scATAC | [GSE201336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201336) | fragments (both matrices, true per-cell coverage) |
| NSCLC scATAC | [GSE270148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270148) | filtered peak matrices (both matrices, true per-cell coverage = total peak counts) |
| HNSC scRNA | [GSE181919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181919) | transcriptional reference |
| TCGA-LIHC | [UCSC Xena](https://xenabrowser.net/datapages/?cohort=TCGA%20Liver%20Cancer%20(LIHC)) | expression + clinical, downstream endpoint |
| TCGA-HNSC | [UCSC Xena](https://xenabrowser.net/datapages/?cohort=TCGA%20Head%20and%20Neck%20Cancer%20(HNSC)) | expression + clinical, downstream endpoint |

Also required: **GENCODE v38** (hg38) annotation GTF, from
<https://www.gencodegenes.org/human/release_38.html>.

## The negative control is constructed here, not deposited

The mislabelled negative control used in the marker read-out is **built by this
repository**, and the error in it is introduced by us. The gene-activity matrix
deposited with GSE227265 carries no row names; we assign its rows by the
coordinate order of a GENCODE annotation, which yields wrong gene identities
because the depositors used a different annotation release. **This is not a defect
in the deposited data.** It is a deliberate corruption, used only to check that the
marker test detects wrong gene identities under every aggregation strategy.

## Expected layout

```
$SCRNA_ATAC_BASE/
├── data/
│   ├── raw/
│   │   ├── TCGA_LIHC/TCGA_LIHC_expression.gz
│   │   ├── TCGA_LIHC/TCGA_LIHC_clinical.tsv
│   │   ├── TCGA_HNSC/TCGA_HNSC_expression.gz
│   │   └── TCGA_HNSC/TCGA_HNSC_clinical.tsv
│   └── processed/
│       ├── integrated/
│       │   ├── atac_hcc_geneactivity.h5ad            # deposited matrix, rows named by coordinate order (negative control)
│       │   ├── atac_hcc_geneactivity_proximal.h5ad   # Script 90
│       │   ├── atac_hcc_geneactivity_enhancer.h5ad   # Script 90
│       │   ├── atac_hcc_final.h5ad                   # cell labels
│       │   └── pseudobulk_reference_v2.csv
│       └── HNSC/pseudobulk_reference_hnsc_v3.csv
└── results/validation/                               # Script 91 writes here
```

## Which test needs what

Each test is deliberately cheap, and none needs the whole tree:

- **Test 1** (`audit_marker_concordance.py`) — the three `atac_hcc_geneactivity*.h5ad`
  files only. Runs in about a minute.
- **Test 2** (`audit_zscale.py`) — the colorectal gene-activity matrix from the
  GSE201336 pipeline.
- **Test 3** (`audit_coad_depth.py`) — the GSE201336 fragment files, tabix-indexed,
  plus the stored per-cell table it is compared against.

The downstream comparison (`scripts/91_...`) additionally needs both TCGA cohorts.

## A caveat on the colorectal cell counts

`scripts/65_coad_atac_full_pipeline.py` admits only the first 10,000 barcodes per
sample (`max_cells=10000`, line 110). The per-cell coverage statistics in Table 1 of
the paper are recomputed from the fragment files and are unaffected, but the cell
counts are properties of that truncated set rather than of the cohort. Anyone
re-running this should raise or remove the cap.
