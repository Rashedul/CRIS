#!/bin/bash 

cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/fastq

# trimmomatic
for f in *1.fastq.gz;
do
fq1=$f
fq2=${f%_1*.fastq.gz*}"_2".fastq.gz
fq=${f%_1*.fastq.gz*}
echo $fq 
ls -l $fq1
ls -l $fq2

pgz_1="trimmomatic/${fq}_1P.fastq.gz"
ugz_1="trimmomatic/${fq}_1U.fastq.gz"
pgz_2="trimmomatic/${fq}_2P.fastq.gz"
ugz_2="trimmomatic/${fq}_2U.fastq.gz"

echo $pgz_1

trimmomatic PE "${fq1}" "${fq1}" "${pgz_1}" "${ugz_1}" "${pgz_2}" "${ugz_2}" \
  SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:adapter_pe.fa:2:40:15

done 


## mixcr for RNA-Seq; --receptor-type IGH try

cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/fastq/trimmomatic


for f in *1P.fastq.gz;
do
fq1=$f
fq2=${f%_1*.fastq.gz*}"_2P".fastq.gz
fq=${f%_1*.fastq.gz*}
echo $fq 
ls -l $fq1
ls -l $fq2

mixcr analyze shotgun \
        --species hs \
        --starting-material rna \
        --only-productive \
        --receptor-type IGH \
        $fq1 $fq2 ./mixcr/$fq

done 



mixcr analyze shotgun \
        --species hs \
        --starting-material rna \
        --only-productive \
        --receptor-type IGH \
        SRR1814049_1.fastq.gz SRR1814049_1.fastq.gz analysis