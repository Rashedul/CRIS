#!/bin/bash 

# check 3 required arguments are provided
if [ -z "$3" ]
  then
    echo "ERROR: Parameters are not supplied in CRIS, please provide a bam file, number of threads and memory."
    exit 1
fi

# check argument $1 is a file
if [ -f "$1" ] 
then
	echo "CRIS is runnig..."
else 
	echo "ERROR: CRIS can't find bam file"
	exit 1
fi

# input parameters
input_bam_file_name="${1}" #bam file name
num_threads=$(echo "${2}" | tr "=" "\t" | awk '{print $2}') #e.g., threads=4
max_memory_assembly=$(echo "${3}" | tr "=" "\t" | awk '{print $2}') #e.g., memory=16G

###############################################################
# Slice bam; PE bam to fastq; de novo assembly of transcripts #
###############################################################

#hg38 ig coordinates
echo "chr14   105600001       106880000
chr14_KI270726v1_random 0       43739
chr15   21710000        22190000
chr16   31950001        33970000
chr16_KI270728v1_random 0       1872759" >Ig_regions.bed

echo "print ig regions"
cat Ig_regions.bed | head

#slice bam by regions of interest. For slice, input file must be coordinate-sorted and indexed
sambamba index $input_bam_file_name
sambamba slice $input_bam_file_name -L Ig_regions.bed -o $input_bam_file_name.slice.bam
sambamba flagstat $input_bam_file_name.slice.bam

picard SamToFastq I=$input_bam_file_name.slice.bam F=$input_bam_file_name.slice.R1.fastq F2=$input_bam_file_name.slice.R2.fastq

printf "Trinity is constructing de novo transcripts...\n\n"

#trinity de novo assembly
ls -l $input_bam_file_name.slice.R1.fastq
ls -l $input_bam_file_name.slice.R2.fastq

Trinity --seqType fq --max_memory $max_memory_assembly --left $input_bam_file_name.slice.R1.fastq --right $input_bam_file_name.slice.R2.fastq --CPU $num_threads --trimmomatic --full_cleanup --no_bowtie;

#######################################################################
# filter IGHV transcripts and identify highly expressed transcripts  #
#######################################################################

printf "\n\n Filtering IGHV transcripts using blastn...\n\n"

# filter IGHV sequnces by blastn 
f=*Trinity.fasta
# get transcript ids having match with ig-genes by blastn
blastn -db ./data/IGHV -query $f -out $input_bam_file_name"_blastn" -num_threads $num_threads -outfmt 6;
# remove duplicated ids
cat $input_bam_file_name"_blastn" | awk '!seen[$1]++' | awk '{print $1}' > $input_bam_file_name"_blastn_HG_IG_remdup.ID"
# filter transcript fasta file by transcript id
cat $f | seqkit grep -f $input_bam_file_name"_blastn_HG_IG_remdup.ID" > $input_bam_file_name"_blastn_HG_IG_remdup.ID.fa"

printf "\n\n running salmon to estimate transcript abundace... \n\n"

# create sailfish index and quantify transcript abundances
fa=$input_bam_file_name"_blastn_HG_IG_remdup.ID.fa"
salmon index -t $fa -i salmon_index -k 25 -p $num_threads
salmon quant -i salmon_index -l "OSR" -1 $input_bam_file_name.slice.R1.fastq -2 $input_bam_file_name.slice.R2.fastq -o salmon_index

# sort transcript ids by tpm values
cat ./salmon_index/quant.sf | sort -nr -k4 | grep -v "Name" | sed  '1i #Name    Length  EffectiveLength TPM     NumReads' >$input_bam_file_name.ig-transcripts.sortedbyTPM.txt
cat $input_bam_file_name.ig-transcripts.sortedbyTPM.txt | awk '{print $1}' | grep -v "Name" >list

# sort fasta file by transcript ids
while read line; do
cat $f | seqkit grep -p $line;
done <list >>$input_bam_file_name.ig-transcripts.sortedbyTPM.fasta

# igblastn for SHM and clonotypes
igblastn -germline_db_V ./data/IGHV -num_alignments_V 3 -germline_db_J ./data/IGHJ -num_alignments_J 3 -germline_db_D ./data/IGHD -num_alignments_D 3 -organism human -query $input_bam_file_name.ig-transcripts.sortedbyTPM.fasta -show_translation -auxiliary_data ./data/human_gl.aux >$input_bam_file_name.IgBLAST_out.txt

rm *slice* *blastn* list *.bam *.bam.bai *Trinity.fasta *bed environment.txt
rm -r salmon_index
rm -r data

cat *ig-transcripts.sortedbyTPM.fasta
cat *ig-transcripts.sortedbyTPM.txt
cat *IgBLAST_out.txt

ls -l *

printf "\n\n #### finished job...\n\n"
