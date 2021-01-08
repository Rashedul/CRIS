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
#wget https://github.com/mozack/vdjer/releases/download/v0.10_reference/vdjer_human_references.tar.gz
#tar -xvf vdjer_human_references.tar

for f in S*Aligned.out.sorted.bam; do 
    echo $f;
    mkdir vdjer.$f;
    vdjer --in $f --rl 90 --t 8 --ins 200 --chain IGH --ref-dir igh > vdjer.sam 2> vdjer.log
    mv vdj_contigs.fa vdjer.$f/vdjer.$f.vdj_contigs.fa;
done

# is ~190 tested for two samples
picard CollectInsertSizeMetrics \
I=SRR1814049Aligned.out.sorted.bam  \
O=insert_size_metrics.txt \
H=insert_size_histogram.pdf \
M=0.5

vdjer --in SRR1814049Aligned.out.sorted.bam --rl 90 --t 8 --ins 190 --chain IGH --ref-dir igh > vdjer.sam 2> vdjer.log

vdjer --in ./vdjer/star.sort.bam --rl 50 --t 8 --ins 175 --chain IGH --ref-dir igh > vdjer.sam 2> vdjer.log



# test run vdjer
git pull https://github.com/mozack/vdjer.git

#test salmon
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/bam/test_run/CRIS/CRIS.SRR1814049_test.bam
salmon index -t SRR1814049_test.bam_blastn_HG_IG_remdup.ID.fa -i salmon_index -k 25 -p 8

salmon quant -i salmon_index -l "OSR" -1 SRR1814049_test.bam.slice.R1.fastq -2 SRR1814049_test.bam.slice.R2.fastq -o transcripts_quant


# cris on star
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/bam_star
for f in *.out.sorted.bam; do 
    echo $f; 
    bash CRIS.sh -inbam $f; 
done

# cris on bwa
cd /projects/epigenomics3/epigenomics3_results/users/rislam/CLL_hg38/IGHV_status_RNAseq/RNAseq17samples_Blachly_PNAS/bam
for f in S*.bam; do 
    echo $f; 
    bash CRIS.sh -inbam $f; 
done


#dbgap
cd /projects/cll_dart/rnaseq/sorted_bwa_bam/

#symlink
while read line; do
    echo $line;
    ln -s /projects/cll_dart/rnaseq/bwa/$line.bam ./;
done </projects/cll_dart/rnaseq/documents/bulk_RNAseq_59runs_completed_ID.tsv


#sort bam by coordinates
for f in *bam;
do echo $f;
sambamba sort -t 40 --tmpdir=/projects/epigenomics3/temp/rislam/ $f;
done

#run cris
while read line; do
    echo $line;
    ls -l $line.sorted.bam;
    bash CRIS.sh -inbam $line.sorted.bam;
done </projects/cll_dart/rnaseq/documents/bulk_RNAseq_59runs_completed_ID.tsv

## igblast output pasring

#get query ids:
less VDJER_html.txt | grep "Query=</b>" | awk '{print $2}' > igBlast_query_ids.txt

#get 8 cols: Top V gene match   Top D gene match    Top J gene match    Chain type  stop codon  V-J frame   Productive  Strand
less VDJER_html.txt | grep -A 1 "Top V gene match" | awk '{gsub("<tr><td>|</td></tr></table>|</td></tr>", ""); print $0}' | awk '{gsub("</td><td>", "\t"); print $0}' | grep "IG" > igBlast_VDJ_8col.txt

#get 5 cols. length      matches     mismatches      gaps    identity(%)
less VDJER_html.txt | grep "<tr><td> Total </td><td>" | awk '{print $6, $8, $10, $12, $14}' > igBlast_blast_5col.txt

#get CDR3 nucleotide and AA seq
less VDJER_html.txt | grep "<tr><td> </td><td>Nucleotide sequence</td><td>Translation</td><td>Start</td><td>End</td><tr><td>CDR3</td><td>" | awk '{gsub("<tr><td> </td><td>Nucleotide sequence</td><td>Translation</td><td>Start</td><td>End</td><tr><td>CDR3</td><td>", ""); print $1, $2}' | awk '{gsub("</td><td>", ""); print $0}' > igBlast_CDR3_nt_AA.txt

#note last 7 transcripts are missing nd and aa info
paste igBlast_query_ids.txt igBlast_VDJ_8col.txt igBlast_blast_5col.txt igBlast_CDR3_nt_AA.txt >igBlast_VDJER_summary.txt








