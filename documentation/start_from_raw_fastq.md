
## Download raw fastq from NeMO and run demultiplex
Landing page:https://assets.nemoarchive.org/dat-aw6czix

FTP site: https://data.nemoarchive.org/bican/grant/BICAN_Mul_PN_Human/salk_ecker/epigenome/nuclei/m3C_seq/human/demultiplexed_fastq/

### Download raw fastq from NeMo
```shell
# download all
sed '1d' BG.NeMO.file_urls.tsv | cut -f 4 > urls.txt
parallel -j 3 --bar wget -c -q -x -nH --cut-dirs=8 -P ./fastq {} < urls.txt
```

```shell
# Download in batch region by region
# CaH,CaB,CaT,Pu,NAC,GP,MGM1,SubTH
region="SubTH"
cd ~/Projects/BICAN
mkdir -p downloaded/${region}
awk -F'\t' -v r="$region" '$3==r {print $4}' ~/Projects/BICAN/BG.NeMO.file_urls.tsv | parallel -j 2 --delay 2 wget -nc -nd -P downloaded/${region} {}
# awk -F'\t' -v r="$region" '$3==r {print "wget -c -P downloaded/"r" "$4}' ~/Projects/BICAN/BG.NeMO.file_urls.tsv
# awk -F'\t' -v r="$region" '$3==r {print $4}' ~/Projects/BICAN/BG.NeMO.file_urls.tsv | xargs -L 1 -I {} wget -nc -nd -P downloaded/${region} {}
# detect bad files:
indir="downloaded/${region}"
echo ${indir}
find "$indir" -name "*.fastq.gz" -print0 \
  | parallel -0 -j 4 'gzip -t {} 2>/dev/null || echo {}' \
  > ${region}.bad_files.txt
# NAC,MGM1 and SubTH are done, no bad file;
```

### Reorganize raw FASTQ files into subfolders for batch demultiplexing
Each plate: 8 lanes * 2 (R1 and R2) = 16 files = 384 cells
Each sample: 16 plates * 16 files = 196 files
```python
import os
import pandas as pd
import glob

region="SubTH"
indir=f"downloaded/{region}"
infiles=glob.glob(os.path.join(indir,"*.fastq.gz"))
print(len(infiles))
df=pd.DataFrame(infiles,columns=['source'])
df.insert(0,"Plate",df.source.apply(lambda x:os.path.basename(x).split('-')[-3]))
df.insert(0,"SampleID",df.Plate.apply(lambda x:'-'.join(x.split('_')[:-2])))
df['dest']=df.apply(lambda x:os.path.join(indir,x.SampleID,os.path.basename(x.source)),axis=1)
for sample in df.SampleID.unique():
    outdir=os.path.join(indir,sample)
    os.makedirs(outdir,exist_ok=True)
for source,dest in df.loc[:,['source','dest']].values:
    os.rename(source,dest)
```

### Run demultiplex
```shell
region="SubTH"
# outdir="mapping"
outdir="single_cell_fastq2"
rm -f run_demultiplex_${region}.sh
for sample in `ls downloaded/${region}`; do
    echo ${sample}
    yap-gcp run_demultiplex --fq_dir="downloaded/${region}/${sample}" --outdir="${outdir}/${region}/${sample}" --gcp=False --n_jobs=8 --print_only True >> run_demultiplex_${region}.sh
    #echo "Pbsgen -submit -c 8 -d 3 -m 150 -name demultiplex_${sample} 'yap-gcp run_demultiplex --fq_dir="downloaded/${region}/${sample}" --outdir="${outdir}/${region}/${sample}" --gcp=False --n_jobs=8'" >> run_demultiplex_${region}.sh
done;
source run_demultiplex_${region}.sh
```

### Run mapping
```shell
yap default-mapping-config --mode m3c --barcode_version V2 --bismark_ref "path_to_your_ref/hg38_ucsc_with_chrL.bismark1" --genome "path_to_your_ref/hg38_ucsc_with_chrL.fa" --chrom_size_path "path_to_your_ref/hg38_ucsc.nochrM.sizes" --hisat3n_dna_ref  "path_to_your_ref/hg38_ucsc_with_chrL" > m3c_config.ini
# choose one between bismark_ref and hisat3n_dna_ref
# see https://github.com/DingWB/cemba_data for detail documentation
yap-gcp run_mapping --workd="mapping/${region}/${sample}" --gcp=False --config_path="m3c_config.ini" --aligner='hisat-3n' --n_jobs=64 --print_only=True
# "{uid}/fastq/{cell_id}-R1.fq.gz" will be detected under workd
```