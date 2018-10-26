#!/bin/csh -f

mkdir fastqs

sed 1d SRR.txt | awk '{print $2}' > SRR.name.txt

foreach i ("`cat SRR.name.txt`")

set ii =  `echo $i | cut -c1-6`

/home/exx/.aspera/connect/bin/ascp -i \
/home/exx/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -QT -l 200m \
 anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/$ii/$i/$i.sra fastqs

end


