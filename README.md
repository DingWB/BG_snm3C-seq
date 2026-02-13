# A Multimodal Single-Cell Epigenomic and 3D Genome Atlas of the Human Basal Ganglia

## 1. Data Availability
### 1.1 Raw fastq files
The raw fastq files (plate-level, before demultiplexing)  are available at:
- NEMO: https://assets.nemoarchive.org/dat-aw6czix; 
- FTP site: https://data.nemoarchive.org/bican/grant/BICAN_Mul_PN_Human/salk_ecker/epigenome/nuclei/m3C_seq/human/demultiplexed_fastq/

### 1.2 cell-level fastq
Cell-level fastq files (after demultiplex) can be downloaded from GEO with accession ID: 

### 1.3 Spatial data
MERFISH transcriptomic data will be available on [BIL](https://doi.org/10.35077/g.1194) soon

### 1.4 bigwig files
bigwig files are available for download on [Basal Ganglia Epigenome Browser](https://basalganglia.epigenomes.net/)

### 1.5 Processed data
Processed data are available on the NCBI Gene Expression Omnibus (GEO) and Figshare
#### GEO
- Pseudobulk allc files (Group level): GEO
- HiC contact files (single-cell level): GEO

#### Figshare (uploading)
- Pseudobulk allc files (Subclass level CG+CH): folder Subclass.allc
- adata: folder adata
- DMG: folder DMG
- DMR: folder DMR
- Enriched motif and TF: folder motif
- Normalized compartment scores across cell types at subclass level: HiC/NormalizedCompartmentScores.tsv
- Diff domain doundary: HiC/diff_boundary.tsv
- Loop at Subclass and Group levels: folder HiC/Subclass.loop and HiC/Group.loop
- Enhancer-Promoter Links: supplementary_tables/TableS7.enhancer_promoter_links.tsv
- TF-Target gene pairs: supplementary_tables/TableS8.Subclass.GRN.xlsx
- Supplementary Tables (for the manuscript): folder supplementary_tables

## 2. Mapping Pipeline
If you prefer to process the data from fastq files, please run dumultiplexing and mapping using our mapping pipeline:
- Local (or Google Cloud) snakemake pipeline: https://github.com/DingWB/cemba_data
- Broad WDL pipeline: https://broadinstitute.github.io/warp/docs/Pipelines/snm3C/README

### 3. Code and jupyter notebooks for clustering, cell type annotations and downstream analysis
- Clustering: See [clustering folder](clustering/)
- Integration between snm3C-seq & sc-RNA: See [integration folder](integration/)
- Call DMR: https://github.com/lhqing/ALLCools
- Motif enrichment: [motif_enrichment folder](motif_enrichment/)
- ABC Model: https://github.com/broadinstitute/ABC-Enhancer-Gene-Prediction/tree/main
- Spatial Transcriptomics Data Analysis: https://github.com/a3klein/BG_snm3c-seq_MERFISH
- LDSC: https://github.com/rosanwang/ldsc-wrapper

### 4. Data exploration & visualization
- (1) Basal Ganglia Epigenome Browser: https://basalganglia.epigenomes.net/
- (2) Human Basal Ganglia (HBG) Anno-J Networked Genome Browser: https://neomorph.salk.edu/hbg/hbg.php

### 5. Citation


## 6. FAQ
### How to normalize the raw methylation fraction?
```python
import os,sys
import pandas as pd
import anndata

adata_path="BG.gene-CHN.h5ad" # example adata path
gene="DRD2" # an example gene
clip_norm_value=10
raw_adata = anndata.read_h5ad(os.path.expanduser(adata_path), backed='r')
adata = raw_adata[:, gene].to_memory()
raw_adata.file.close()
cols = adata.obs.columns.tolist()
if 'prior_mean' in cols:
    na_sum = adata.to_df().isna().sum().sum()
    adata.X = adata.X / adata.obs.prior_mean.values[:, None]  # range = [0,1,10]
    if not clip_norm_value is None:
        if issparse(adata.X):
            X=adata.X.toarray()
        else:
            X=adata.X
        adata.X = np.clip(X, None, clip_norm_value)
    adata.uns['normalize_per_cell'] = True
adata
```

