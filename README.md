# Multimodal single-cell epigenomic and 3D genome atlas of the human basal ganglia

## 1. Data Availability
For instructions on downloading and processing the raw FASTQ files, see [start_from_raw_fastq.md](documentation/start_from_raw_fastq.md).

### 1.1 Raw FASTQ files
The raw FASTQ files (plate-level, before demultiplexing) are available at:
- NeMO: https://assets.nemoarchive.org/dat-aw6czix
- FTP site: https://data.nemoarchive.org/bican/grant/BICAN_Mul_PN_Human/salk_ecker/epigenome/nuclei/m3C_seq/human/demultiplexed_fastq/

### 1.2 Cell-level FASTQ files
Cell-level FASTQ files (after demultiplexing) can be downloaded from GEO under accession [GSE320293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320293) or from [AWS Open Data](https://registry.opendata.aws/salk-human-brain-snm3c-seq).

### 1.3 Spatial data
MERFISH transcriptomic data is available on [BIL](https://doi.org/10.35077/g.1194)

### 1.4 BigWig, Hi-C, and track files
BigWig, Hi-C, and other track files can be visualized on the [Basal Ganglia Epigenome Browser](https://basalganglia.epigenomes.net/) or downloaded from the [FTP site](https://epigenome.wustl.edu/basal-ganglia-epigenome/tracks/eckerlab/).

### 1.5 Processed data
Processed data are available through the NCBI Gene Expression Omnibus (GEO) and [AWS Open Data](https://registry.opendata.aws/salk-human-brain-snm3c-seq).
#### GEO: [GSE320293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320293)
- Pseudobulk allc files (Group level)
- HiC contact files (single-cell level)
- Single-cell level fastq files (uploading)

#### AWS Open Data
- single_cell_fastq: raw single-cell FASTQ files
- WGS: raw whole-genome sequencing FASTQ files
- allc: single-cell and pseudobulk allc.tsv.gz files
- hic: pseudobulk .hic files
- mcool: pseudobulk 3D genome contacts in mcool format
- loop.Q.cool: imputed 3D genome contact files in cool format at 10-kb resolution
- h5ad: AnnData files containing gene expression, 100-kb mCG and mCH, and 100- and 25-kb 3D genome embeddings
- ABC: enhancer-promoter links predicted by the ABC model

## 2. Mapping Pipeline
To process the FASTQ files, run demultiplexing and mapping using one of these pipelines:
- Local (or Google Cloud) snakemake pipeline: https://github.com/DingWB/cemba_data
- Broad WDL pipeline: https://broadinstitute.github.io/warp/docs/Pipelines/snm3C/README

For an example workflow that downloads raw FASTQ files and runs mapping, see [start_from_raw_fastq.md](documentation/start_from_raw_fastq.md). 

## 3. Code and Jupyter Notebooks for Clustering, Cell-Type Annotation, and Downstream Analysis
- Clustering: See [clustering folder](clustering/)
- Integration between snm3C-seq & sc-RNA: See [integration folder](integration/)
- Call DMR: https://github.com/lhqing/ALLCools
- HiC: See [HiC folder](HiC/)
- Motif enrichment: [motif_enrichment folder](motif_enrichment/)
- ABC Model: [ABC folder](ABC/)
- Spatial Transcriptomics Data Analysis: https://github.com/a3klein/BG_snm3c-seq_MERFISH
- LDSC: https://github.com/rosanwang/ldsc-wrapper

## 4. Data Exploration and Visualization
- (1) Basal Ganglia Epigenome Browser: https://basalganglia.epigenomes.net/
- (2) Human Basal Ganglia (HBG) Anno-J Networked Genome Browser: https://neomorph.salk.edu/hbg/hbg.php
- (3) Single-cell UMAP, boxplot, dot heatmap and Multiomic track view: https://neomorph.salk.edu/SCMDAP/BasalGanglia

## 5. Citation
Ding, W., Klein, A., Baez-Becerra, C. T., Rink, J. A., Bartlett, A., Zeng, Q., ... & Ecker, J. R. (2026). [A Multimodal Single-Cell Epigenomic and 3D Genome Atlas of the Human Basal Ganglia](https://www.biorxiv.org/content/10.64898/2026.02.12.705594v1). bioRxiv, 2026-02.

## 6. FAQ
### How to normalize the raw methylation fraction?
```python
import os
import numpy as np
import pandas as pd
import anndata
from scipy.sparse import issparse

adata_path = "BG.gene-CHN.h5ad"  # example AnnData path
gene = "DRD2"  # example gene
clip_norm_value = 10
raw_adata = anndata.read_h5ad(os.path.expanduser(adata_path), backed='r')
adata = raw_adata[:, gene].to_memory()
raw_adata.file.close()
cols = adata.obs.columns.tolist()
if 'prior_mean' in cols:
    adata.X = adata.X / adata.obs.prior_mean.values[:, None]  # range = [0,1,10]
    if clip_norm_value is not None:
        if issparse(adata.X):
            X = adata.X.toarray()
        else:
            X = adata.X
        adata.X = np.clip(X, None, clip_norm_value)
    adata.uns['normalize_per_cell'] = True
adata
```

### cell id in metadata does not match the cell id from the mapping pipeline
The Broad mapping pipeline generates cell IDs slightly differently. Please use the following code to convert the cell IDs from the Broad pipeline into new_cell_id:
```python
import pandas as pd
def _parse_cell_id_v2(cell_id):
	plate, multiplex_group, pcr_index, random_index = cell_id.split('-')
	# 384 pos
	col384 = int(random_index[1:]) - 1
	row384 = ord(random_index[0]) - 65  # convert A-P to 0-15
	resolved_multiplex_group=(col384 % 12) // 2 + 1 if 'unknow' not in random_index.lower() else 'NA'
	uid=plate+'-'+str(resolved_multiplex_group)+'-'+pcr_index
	new_cell_id = uid + '-' + random_index
	record = pd.Series({
		'new_cell_id':new_cell_id,
		'Plate': plate,
		'PCRIndex': pcr_index,
		'MultiplexGroup': multiplex_group,
		'ResolvedMultiplexGroup': resolved_multiplex_group,
		'RandomIndex': random_index,
		'Col384': col384,
		'Row384': row384
	})
	return record

metadata=pd.read_csv("annotations.tsv",sep='\t',index_col=0)
cell_ids=metadata.index.tolist()
plate_info = pd.DataFrame([_parse_cell_id_v2(cell_id) for cell_id in cell_ids],index=cell_ids)
# use new_cell_id
```

## 7. Citations
If you use the snm3C-seq or MERFISH datasets generated by this study, or SCMDAP, please cite:

- Ding W, Klein A, Báez-Becerra C ...
[Multimodal single-cell epigenomic and 3D genome atlas of the human basal ganglia](https://www.cell.com/cell/fulltext/S0092-8674(26)01269-9). Cell, 189, 7386-7405.e15 </br>
PMID: 43020023</br>
DOI: 10.1016/j.cell.2026.10.019

If you use the Basal Ganglia Epigenome Portal (https://basalganglia.epigenomes.net), please <b>**also**</b> cite: 

- Zhang W, Ding W, Li K ...
[An integrated single-cell and epigenomic resource for comparative analysis of the basal ganglia from BICAN](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(26)00268-7). Cell Reports Methods, 2026; 0 </br>
PMID: 43020065 </br>
DOI: 10.1016/j.crmeth.2026.101567

If you use data from SCMDAP or the Basal Ganglia Epigenome Portal, please also cite the relevant [BICAN publication](https://www.cell.com/consortium/bican-basal-ganglia-atlas).
