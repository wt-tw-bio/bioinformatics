wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
#下载完成之后解压
gunzip ko_list.gz
tar -xzvf profiles.tar.gz 
wget https://www.genome.jp/ftp/tools/kofam_scan/kofam_scan-1.3.0.tar.gz
tar -xzvf kofam_scan-1.3.0.tar.gz
echo export PATH=/home/wtao/projects/bacteria/kofam_scan-1.3.0:\$PATH >> ~/.bashrc
source ~/.bashrc
conda activate kegg
conda install -c conda-forge ruby
conda install -c bioconda hmmer
conda install -c conda-forge parallel
cd kofam_scan-1.3.0
cp config-template.yml config.yml
exec_annotation -o CKP28_keggid.txt --cpu 20 --format mapper -E 1e-5 /home/wtao/projects/czh/P24052701/bacteria/outdir/Prokka/CKP28/CKP28.faa