# 基因组比对
conda install -c bioconda mummer
nucmer --prefix=genome_comparison ./参考基因组/ncbi_dataset_k2044/ncbi_dataset/data/GCA_000009885.1/GCA_000009885.1_ASM988v1_genomic.fna ./P24052701/bacteria/outdir/Prokka/CKP28/CKP28.fna
#  --prefix=genome_comparison: 设置输出文件的前缀
#  第一个文件路径是参考基因组
#  第二个文件路径是要比对的目标基因组
#  输出将生成 .delta 文件，包含比对结果
show-coords -q -c -l genome_comparison.delta > genome_comparison.coords
#  将 delta 文件转换为易读的坐标格式：
#  -q: 基于查询序列(query)的排序方式展示比对结果
#  -c: 显示覆盖度信息
#  -l: 显示序列长度
show-diff -q genome_comparison.delta > genome_differences3.txt
#  分析并输出基因组间的差异：
#  -q: 基于查询序列的分析模式
#  输出结果包含插入、缺失和重排等结构变异信息

# 毒力因子预测
#  下载blast+
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.16.0+-x64-linux.tar.gz 
mv ncbi-blast-2.16.0+  blast
#  export PATH=/work/run/projects/wtao/czh/blast/bin:$PATH
source ~/.bashrc
#  下载毒力因子数据库：https://www.mgc.ac.cn/VFs/download.htm
gzip -d VFDB_setB_pro.fas.gz
#  构建数据库
makeblastdb -in VFDB_setB_pro.fas -dbtype prot -out vfdb
#  与数据库比对
blastp -query P24052701/bacteria/outdir/Prokka/CKP28/CKP28.faa  -db vfdb -out vfdb_results.txt -evalue 1e-5 -outfmt "6 qseqid sseqid qstart qend length qlen slen evalue bitscore stitle" -num_threads 20

# COG注释
#  diamond跟blast+差不多，速度更快
#  diamond下载：
wget -c https://github.com/bbuchfink/diamond/releases/download/v2.1.11/diamond-linux64.tar.gz
tar -zxvf diamond-linux64.tar.gz 
#  COG数据库下载
wget -c https://ftp.ncbi.nlm.nih.gov/pub/COG/COG2024/data/COGorg24.faa.gz
gunzip COGorg24.faa.gz
#  run
./diamond help
../diamond/diamond makedb --in COGorg24.faa -d cog_db
../diamond/diamond blastp -d cog_db -q ../P24052701/bacteria/outdir/Prokka/CKP28/CKP28.faa -out cog_results.txt -e 1e-5 -k 1 -p 20 --outfmt 6 qseqid sseqid qstart qend length qlen slen evalue bitscore stitle
../diamond/diamond blastp \
-d cog_db \
-q ../P24052701/bacteria/outdir/Prokka/CKP28/CKP28.faa \
-o cog_results.txt \
-e 1e-5 \
-k 1 \
-p 20

# circos
#  下载
conda install circos
circos -v
#  circos文件地址
conda list circos | grep circos
#  run
circos -conf circos.cnof