# Snakemake-CHIPseq

## Installation related software in CentOS 7

#### install MACS1 and MACS2 through python-pip

curl "https://bootstrap.pypa.io/get-pip.py" -o "get-pip.py"

python get-pip.py

pip -V

pip install MACS2

pip install MACS=1.4.2

#### download miniconda python3 64bit and install following software
conda config --add channels defaults

conda config --add channels conda-forge

conda config --add channels bioconda

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

install.packages("./spp_1.14.tar.gz")   #fialed with error: cannot find Boost headers version >= 1.41.0

sudo apt-get install libboost1.54-all-dev  #yum install boost-devel -y

R

install.packages("./spp_1.14.tar.gz")

still failed in ubuntu 18.04

### https://ubuntuforums.org/showthread.php?t=2363177
sudo apt-get --purge remove libboost-dev libboost-doc
sudo apt-get --purge remove libboost-dev
sudo apt-get --purge remove libboost-all-dev
sudo apt autoremove

wget https://dl.bintray.com/boostorg/release/1.68.0/source/boost_1_68_0.tar.bz2

tar xvjf  boost_1_68_0.tar.bz2

cd boost_1_68_0

./bootstrap.sh --prefix=/usr/local

sudo apt-get update
sudo apt-get install build-essential g++ python-dev autotools-dev libicu-dev build-essential libbz2-dev libboost-all-dev

./bootstrap.sh --prefix=/usr/local

./b2

sudo ./b2 install

dpkg -s libboost-dev | grep 'Version'

#### prepare ROSE
wget https://bitbucket.org/young_computation/rose/get/1a9bb86b5464.zip

gunzip 1a9bb86b5464.zip


#### install parallel

wget -O - pi.dk/3 # (wget -O - pi.dk/3 || curl pi.dk/3/ || fetch -o - pi.dk/3)  

bash


### Installation related software in Ubuntu 18.04
wget https://github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz

tar -zxvf MACS-1.4.2-1.tar.gz

cd MACS-1.4.2-1

python setup.py install

macs14 -help


wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

sh Miniconda3-latest-Linux-x86_64.sh

bash 

conda config --add channels defaults

conda config --add channels conda-forge

conda config --add channels bioconda

conda add -c bioconda MACS2   #macs2 --help

conda install multiqc

conda install fastqc

conda install bowtie

conda install samblaster

conda install samtools

conda install deeptools

conda install R



conda install -c bioconda sambamba

conda install -c bioconda sratoolkit   ##if failed, soinstall manually

wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.4.1/sratoolkit.2.4.1-ubuntu64.tar.gz

tar xzvf sratoolkit.2.4.1-ubuntu64.tar.gz

export PATH=$PATH:/home/qcai1/Downloads/sratoolkit.2.4.1-ubuntu64/bin 

conda install -c bioconda -c conda-forge snakemake  #failed

Solving environment: failed

UnsatisfiableError: The following specifications were found to be in conflict:
  - snakemake
  - subprocess32
Use "conda info <package>" to see the dependencies for each package.

sudo apt install python3-pip

pip3 install snakemake  #failed

sudo apt install snakemake  #works

### Mapping in HPC without snakemake

### reference genome file
ftp://ftp.ensembl.org/pub/release-81/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
ftp://ftp.ensembl.org/pub/release-81/gtf/homo_sapiens/Homo_sapiens.GRCh38.81.gtf.gz
bowtie-build /home/exx/Documents/pyflow-ChIPseq-master/ht38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa ht38


#### download prebluit indexs 
http://bowtie-bio.sourceforge.net/tutorial.shtml


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
  cd fastqs

  #### convert sra to fastqs and compress to .gz files

  #### you can use a for loop to fastq-dump the downloaded sra files.
find *sra| parallel -j 4  fastq-dump {}

find *fastq | parallel -j 4  bgzip {}

  #### save some space
rm *sra

  #### sample.list making
python sample2json_C.py --fastq_dir fastqs/ --meta SRR.txt


### rerun
snakemake 04aln_downsample/MOLM-14_DMSO1_5_BRD4-downsample.sorted.bam



http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html

nakemake --snakefile Snakefile -np --dag | dot -T png >  t.png
nakemake --snakefile Snakefile -np --rulegraph | dot -T png >  t2.png



