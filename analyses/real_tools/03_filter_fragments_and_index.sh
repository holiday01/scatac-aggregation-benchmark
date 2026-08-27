#!/bin/bash
# Subset the GSE227265 all-sample fragment file (1.5e9 lines, all unfiltered barcodes of 13 samples) to the
# 12,029 annotated HCC cells. This is the equivalent of Signac::FilterCells() (documented for exactly this
# purpose) done with grep/mawk because FilterCells reads the whole file in single-threaded R.
# Then bgzip + tabix (pysam) so that Signac's CreateFragmentObject/FeatureMatrix can use the index.
set -euo pipefail
SCR=/tmp/claude-1000/-mnt-10t-scrna-atac/b849b398-690b-4bff-b242-236b530a2a05/scratchpad
SRC="${GSE227265_FRAGMENTS:-data/raw/GSE227265/GSE227265_fragments_AllSamples.tsv}"
cd $SCR/hcc_frag
echo "$(date +%T) filtering"
LC_ALL=C grep -F -w -f barcodes.txt $SRC | mawk 'BEGIN{while((getline l < "barcodes.txt")>0) s[l]=1} ($4 in s)' > hcc_12029.fragments.tsv
echo "$(date +%T) lines: $(wc -l < hcc_12029.fragments.tsv)"
python3 - <<'PY'
import pysam
pysam.tabix_compress("hcc_12029.fragments.tsv", "hcc_12029.fragments.tsv.gz", force=True)
pysam.tabix_index("hcc_12029.fragments.tsv.gz", preset="bed", force=True)
print("bgzipped + indexed")
PY
rm -f hcc_12029.fragments.tsv
awk -F'\t' '$3=="gene"' "${GENCODE_GTF:-data/reference/gencode.v38.annotation.gtf}" > gencode.v38.genes.gtf
echo "$(date +%T) gene records: $(wc -l < gencode.v38.genes.gtf)"
ls -la
