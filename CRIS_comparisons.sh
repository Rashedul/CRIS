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

#test
mixcr analyze shotgun \
        --species hs \
        --starting-material rna \
        --only-productive \
        --receptor-type IGH \
        SRR1814049_1.fastq.gz SRR1814049_1.fastq.gz analysis

## star alignment

#convert gz to fastq 
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/fastq/
cd /projects/epigenomics3/temp/rislam/Blachly_fastq/
for f in *gz; do 
    echo $f;
    gunzip $f;
done

cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/bam_star

# make star index
mat=/projects/epigenomics3/epigenomics3_results/users/rislam/tool/rMATS.3.2.5
ens=/projects/blueprint_CLL_dart2/analysis/rnaSeq/rMATS/ens/
hg=/projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/genome_feature/hg38_no_alt.fa
cd $ens
$py2/bin/STAR --runThreadN 40 --runMode genomeGenerate --genomeDir $ens/star_index_GRCh38.82 --genomeFastaFiles $hg --sjdbGTFfile $mat/gtf/Homo_sapiens.Ensembl.GRCh38.82.gtf

##alignment 

fq=/projects/epigenomics3/temp/rislam/Blachly_fastq/
ens=/projects/blueprint_CLL_dart2/analysis/rnaSeq/rMATS/ens/
star=/projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/bam_star/
cd $star

for f in  $fq/*_1.fastq;
do 
BamName=$(basename $f _1.fastq);
echo $BamName;
ls -l $fq/$BamName"_1.fastq"
ls -l $fq/$BamName"_2.fastq"
STAR --runMode alignReads --outSAMtype BAM Unsorted --outSAMunmapped Within --outFileNamePrefix $BamName --runThreadN 40 --genomeDir $ens/star_index_GRCh38.82 --readFilesIn $fq/$BamName"_1.fastq" $fq/$BamName"_2.fastq";
done


#sort bam by coordinates
for f in *bam;
do echo $f;
sambamba sort -t 40 --tmpdir=$star $f;
done


#vdjer
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/bam_star/
wget https://github.com/mozack/vdjer/releases/download/v0.10_reference/vdjer_human_references.tar.gz
tar -xvf vdjer_human_references.tar

#vdjer --in star.sort.bam --t 8 --ins 175 --chain IGH --ref-dir vdjer_human_references/igh > vdjer.sam
vdjer --in SRR1814049Aligned.out.sorted.bam --t 8 --ins 175 --chain IGH --ref-dir igh > vdjer.sam


#test salmon
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/bam/test_run/CRIS/CRIS.SRR1814049_test.bam
salmon index -t SRR1814049_test.bam_blastn_HG_IG_remdup.ID.fa -i salmon_index -k 25 -p 8

salmon quant -i salmon_index -l "OSR" -1 SRR1814049_test.bam.slice.R1.fastq -2 SRR1814049_test.bam.slice.R2.fastq -o transcripts_quant














