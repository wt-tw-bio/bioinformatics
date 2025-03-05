# 创建环境
conda create bacteria
# 下载 nextflow
conda install nextflow
# 下载 nf-core
conda install nf-core
nf-core download nf-core/bacass
tar -zxvf nf-core-bacass_2.4.0.tar.gz
#run
nextflow run ./nf-core-bacass_2.4.0/2_4_0 \
 --input samplesheet.csv \
 -profile singularity \
 --outdir outdir \
 --kraken2db 'https://genome-idx.s3.amazonaws.com/kraken/k2_standard_8gb_20210517.tar.gz' \
 --kmerfinderdb /home/wtao/projects/czh/kmerfinder_db/database/bacteria/ \
 --ncbi_assembly_metadata assembly_summary.txt \
 -resume



