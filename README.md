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

save following scripts as download_SRR.csh (adjust ascp address)

##############

#!/bin/csh -f

mkdir fastqs

sed 1d SRR.txt | awk '{print $2}' > SRR.name.txt

foreach i ("`cat SRR.name.txt`")

    set ii =  `echo $i | cut -c1-6`
    
    /home/exx/.aspera/connect/bin/ascp -i \
    
        /home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m \
        
        anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/$ii/$i/$i.sra fastqs
end

 #########################
 
  tcsh download_SRR.csh
 
 # all the sra files will be downloaded in the current fastqs folder.
 convert sra to fastqs and compress to .gz files

## you can use a for loop to fastq-dump the downloaded sra files.
find *sra| parallel -j 4  fastq-dump {}

find *fastq | parallel -j 4  bgzip {}

## save some space
rm *sra

## sample.list making
python sample2json_C.py --fastq_dir fastqs/ --meta SRR.txt


