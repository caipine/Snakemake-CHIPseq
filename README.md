# Snakemake-CHIPseq

This is a pipeline for CHIPseq analysis which was originally built by crazyhottommy (https://github.com/crazyhottommy/pyflow-ChIPseq), here I will do some modifcation to fit my  work statation.

## Installation related software in CentOS 7

#### install MACS1 and MACS2 through python-pip

curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"

python get-pip.py

pip -V

pip install MACS2

pip install MACS=1.4.2

#### download miniconda python3 64bit and install following software
conda config --add channels defaults

conda config --add channels bioconda

conda config --add channels conda-forge

conda install multiqc

conda install fastqc

conda install bowtie

conda install samblaster

conda install samtools

conda install deeptools

conda install R

conda install -c bioconda -c conda-forge snakemake

conda install sratoolkit

conda install -c bioconda sambamba


#### install phantompeakqualtools in R

git clone https://github.com/kundajelab/phantompeakqualtools

R

install.packages("snow", repos="http://cran.us.r-project.org")

install.packages("snowfall", repos="http://cran.us.r-project.org")

install.packages("bitops", repos="http://cran.us.r-project.org")

install.packages("caTools", repos="http://cran.us.r-project.org")

source("http://bioconductor.org/biocLite.R")

biocLite("Rsamtools",suppressUpdates=TRUE)

install.packages("./spp_1.14.tar.gz")

#### prepare ROSE


#### install parallel

wget -O - pi.dk/3 # (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - pi.dk/3)  

bash


### Installation related software in Ubuntu 18.04

### Mapping in HPC without snakemake

### reference genome file
ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
ftp://ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz
bowtie-build /home/exx/Documents/pyflow-ChIPseq-master/ht38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa ht38

### Download SRA data
save follow text as SRR.txt

sample_name fastq_name  factor
MOLM-14_DMSO1_5 SRR2518123   BRD4
MOLM-14_DMSO1_5 SRR2518124  Input
MOLM-14_DMSO2_6 SRR2518125  BRD4
MOLM-14_DMSO2_6 SRR2518126  Input

/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518123/SRR2518123.sra ./

/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518124/SRR2518124.sra ./

/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518125/SRR2518125.sra ./

/home/exx/.aspera/connect/bin/ascp -i  /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR251/SRR2518126/SRR2518126.sra ./

### sample.list making



