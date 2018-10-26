# Snakemake-CHIPseq

a pipeline for CHIPseq analysis

(this pipeline was originally  built by crazyhottommy, here I will add some messager about using experience.)

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

### Download SRA data

### sample.list making



