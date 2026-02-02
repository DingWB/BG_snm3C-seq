## Tool to run motif enrichment on hypo- or hyper-DMRs:
https://github.com/aertslab/pycistarget

## Step1: Scan the DMR sequence and create database
```shell
# install softwares
cd ~/Software/ && git clone https://github.com/aertslab/create_cisTarget_databases.git
cd ~/Software/ && git clone https://github.com/weng-lab/cluster-buster && cd cluster-buster && make cbust
chmod 777 ~/Software/cluster-buster/cbust
pip install flatbuffers

# Creating cisTarget databases: 
# https://pycistarget.readthedocs.io/en/latest/pycistarget_chip-RTD.html#A.-Creating-your-DEM-databases
genome_fasta="${HOME}/Ref/hg38/hg38_ucsc.fa"
path_to_motif_collection="${HOME}/Ref/hg38/annotations/v10nr_clust_public/singletons" #https://resources.aertslab.org/cistarget/motif_collections/
motif_list="${HOME}/Ref/hg38/annotations/v10nr_clust_public/motif_filenames.tsv"
create_cistarget_databases_dir="${HOME}/Software/create_cisTarget_databases"
cbust_path="${HOME}/Software/cluster-buster/cbust"
n_cpu=128 # Srun 128 900 highmem 1
rm -f create_database.sh

folder="Other"
for group in `ls ${folder}/`; do
    echo ${group}
    region_bed="${folder}/${group}/DMR.bed" # your DMR bed file, could be DMR slop 250 bp or the original DMR region
    region_fasta="${folder}/${group}/DMR.fa" # extracted DMR sequence
    output_prefix="${folder}/${group}/DMR.dem_db"
    if ! [ -f "${region_fasta}" ]; then
        bedtools getfasta -fi ${genome_fasta} -bed ${region_bed} > ${region_fasta}
    fi;
    #### Score the motifs: generate motifs_vs_regions.scores.feather, regions_vs_motifs.scores.feather and regions_vs_motifs.rankings.feather
    if ! [ -f "${folder}/${group}/DMR.dem_db.regions_vs_motifs.rankings.feather" ]; then
        echo "${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
-f ${region_fasta} -M ${path_to_motif_collection} -c ${cbust_path} \
-m ${motif_list} -o ${output_prefix} -t ${n_cpu} -l -s 0" >> create_database.sh
    fi;
    ### Create rankings for cisTarget database only
    # motifs_vs_regions_scores_feather="${folder}/${group}/AllDMR.dem_db.motifs_vs_regions.scores.feather"
    # ${create_cistarget_databases_dir}/convert_motifs_or_tracks_vs_regions_or_genes_scores_to_rankings_cistarget_dbs.py -i ${motifs_vs_regions_scores_feather} -s 0
done;

source create_database.sh
```

## Step2: Run pycistarget
```shell
# cistarget_db="${HOME}/Ref/hg38/annotations/cistarget_database/hg38_screen_v10_clust.regions_vs_motifs.rankings.feather"
# dem_db="${HOME}/Ref/hg38/annotations/cistarget_database/hg38_screen_v10_clust.regions_vs_motifs.scores.feather"
# database of motif scores for genomic regions: regions_vs_motifs.scores.feather
path_to_motif_annotations="${HOME}/Ref/hg38/annotations/cistarget_database/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
folder="Other"
rm -f ${folder}.hypo_cistarget_cmd ${folder}.hypo_dem_cmd ${folder}.hyper_cistarget_cmd ${folder}.hyper_dem_cmd
# cell_type_specific_bed or bed
for group in `ls ${folder}/`; do
    rm -rf ${folder}/${group}/*/hypo ${folder}/${group}/*/hyper
    # rm -rf ${folder}/${group}/*/hypo_dem ${folder}/${group}/*/hyper_dem
    cistarget_db="${folder}/${group}/DMR.dem_db.regions_vs_motifs.rankings.feather" # pd.read_feather, columns are: region1, region2,..., motifs
    dem_db="${folder}/${group}/DMR.dem_db.regions_vs_motifs.scores.feather"
    if [ -d "${folder}/${group}/bed" ]; then
        echo ${group}
        for hypo_bedfile in `ls ${folder}/${group}/bed/*.hypo.dmr.bed`; do
            hyper_bedfile=${hypo_bedfile/.hypo.dmr./.hyper.dmr.}
            hypo_background=$(echo $(ls ${folder}/${group}/bed/*.hypo.dmr.bed | grep -v "$(basename ${hypo_bedfile})"))
            hyper_background=$(echo $(ls ${folder}/${group}/bed/*.hyper.dmr.bed | grep -v "$(basename ${hyper_bedfile})"))
            subgroup=$(echo ${hypo_bedfile/.hypo.dmr.bed/} | cut -d "/" -f 4)
            echo ${hypo_bedfile} ${hyper_bedfile}
            mkdir -p ${folder}/${group}/${subgroup}/hypo ${folder}/${group}/${subgroup}/hyper
            n_hypo=$(wc -l ${hypo_bedfile} | cut -d " " -f 1)
            n_hyper=$(wc -l ${hyper_bedfile} | cut -d " " -f 1)
            n_hypo_background=$(( n_hypo * 10 ))
            n_hyper_background=$(( n_hyper * 10 ))
            echo ${n_hypo} ${n_hyper} ${n_hypo_background} ${n_hyper_background}
            echo "pycistarget cistarget --cistarget_db_fname ${cistarget_db} --bed_fname ${hypo_bedfile} --output_folder ${folder}/${group}/${subgroup}/hypo --specie homo_sapiens --path_to_motif_annotations ${path_to_motif_annotations} --output_mode hdf5 --write_html" >> ${folder}.hypo_cistarget_cmd
            echo "pycistarget cistarget --cistarget_db_fname ${cistarget_db} --bed_fname ${hyper_bedfile} --output_folder ${folder}/${group}/${subgroup}/hyper --specie homo_sapiens --path_to_motif_annotations ${path_to_motif_annotations} --output_mode hdf5 --write_html" >> ${folder}.hyper_cistarget_cmd
        done;
    fi;
done;
cat ${folder}.hypo_cistarget_cmd | parallel -j 4
cat ${folder}.hyper_cistarget_cmd | parallel  -j 4

ls ${folder}/*/*/hypo/motif_enrichment*hypo.dmr.html | wc -l
ls ${folder}/*/*/hyper/motif_enrichment*hyper.dmr.html | wc -l
ls ${folder}/*/*/hypo_dem/motif_enrichment_dem_*.html | wc -l
ls ${folder}/*/*/hyper_dem/motif_enrichment_dem_*.html | wc -l
```