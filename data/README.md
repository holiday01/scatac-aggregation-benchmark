# Data

No data files are stored in this repository. Everything analysed is public.

| accession | files used |
|---|---|
| [GSE227265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE227265) | fragment files (13 hepatocellular tumours); the deposited gene-activity matrix is used only to build the deliberately mislabelled negative control |
| [GSE201336](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201336) | fragment files for GSM6058771, GSM6058772, GSM6058774, GSM6058775, GSM6058779, GSM6058785 |
| [GSE270148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270148) | filtered peak-by-cell matrices, GSM8069376–GSM8069384 |
| 10x Genomics `pbmc_granulocyte_sorted_10k` | `filtered_feature_bc_matrix.h5` and `atac_fragments.tsv.gz` (+ index) |

Scripts resolve paths from `$SCRNA_ATAC_BASE` (default `/mnt/10t/scrna_atac`).

## Verified checksums, GSE270148 peak matrices

Re-downloaded from GEO on 2026-08-27 and compared with the working copies. Eight
files were byte-identical; `GSM8069383` was not, and the working copy was
replaced (see `../CHANGELOG.md`).

| sample | patient | bytes | MD5 |
|---|---|---|---|
| GSM8069376 | Lu881 | 42,702,046 | 8d51b78d78a628209859e52c05b5ec8d |
| GSM8069377 | Lu883 | 4,359,847 | ae5d923a3e87a8329b1b1bdca5d9b399 |
| GSM8069378 | Lu885 | 81,758,220 | 53458fecffb73e059eb3dc5049952385 |
| GSM8069379 | Lu898 | 40,101,985 | f4dbb2bb5984c73c284e67a422183387 |
| GSM8069380 | Lu912 | 89,096,964 | c582ba09d00e9c00de0c1196a5e2a56f |
| GSM8069381 | Lu929 | 23,367,057 | 6bb23b324a37c541b10366ae6faa3d0e |
| GSM8069382 | Lu931 | 31,338,664 | 861275f5e56dafbe5cae4687eca0fba3 |
| GSM8069383 | Lu934 | 36,768,155 | 826ccbbc53226c88a6469e17a2a2e6ff |
| GSM8069384 | Lu936 | 77,342,602 | 261836bc60ca46d797e995b2dfb242c4 |
