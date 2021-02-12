#!/bin/bash 
if [ "$1" == "-h" ] ; then
        echo -e "CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq.

Usage: `basename $0` -inbam <input_bam_file> -outdir <output_directory> -threads <num_threads> -memory <max_memory_assembly>
                        <input_bam_file>: (required) bam file must be aligned to hg38 genome build, coordinate-sorted and indexed
                        <output_directory>: (optional) full path of output directory or output files will be written in current directory
                        <num_threads>: (optional) number of threads; default 4
                        <max_memory_assembly>: (optional) maximum memory in G (gigabyte) allowed for assembly; default 4G"
exit 0
fi

#set defaults
num_threads=4
max_memory_assembly=4G
output_directory=../

while [ $# -gt 0 ]
do
        case "$1" in
        -inbam) input_bam_file="$2"; shift;;
        -outdir) output_directory="$2"; shift;;
        -threads) num_threads="$2"; shift;;
        -memory) max_memory_assembly="$2"; shift;;
        esac
        shift
done 

###############################################################
# Slice bam; PE bam to fastq; de novo assembly of transcripts #
###############################################################
input_bam_file_name=$(basename $input_bam_file)
mkdir -p CRIS_out.$input_bam_file_name

#hg38 ig coordinates
echo "chr14   105600001       106880000
chr14_KI270726v1_random 0       43739
chr15   21710000        22190000
chr16   31950001        33970000
chr16_KI270728v1_random 0       1872759" >CRIS_out.$input_bam_file_name/Ig_regions.bed

#slice bam by regions of interest. For slice, input file must be coordinate-sorted and indexed
sambamba slice $input_bam_file -L CRIS_out.$input_bam_file_name/Ig_regions.bed -o CRIS_out.$input_bam_file_name/$input_bam_file_name.slice.bam

cd CRIS_out.$input_bam_file_name
picard SamToFastq I=$input_bam_file_name.slice.bam F=$input_bam_file_name.slice.R1.fastq F2=$input_bam_file_name.slice.R2.fastq

printf "Trinity is constructing de novo transcripts...\n\n"

#trinity de novo assembly
Trinity --seqType fq --max_memory $max_memory_assembly --left $input_bam_file_name.slice.R1.fastq --right $input_bam_file_name.slice.R2.fastq --CPU $num_threads --trimmomatic --full_cleanup --no_bowtie;

#######################################################################
#                  make blast db for igh sequences                    #
#######################################################################

# # download germline IGHV sequences from IMGT 
# cd ../data
# curl -O http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta
# curl -O http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta
# curl -O http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHJ.fasta

# # modify IMGT fasta for blast compatibility
# edit_imgt_file.pl IGHV.fasta > formatted_IGHV.fasta
# edit_imgt_file.pl IGHD.fasta > formatted_IGHD.fasta
# edit_imgt_file.pl IGHJ.fasta > formatted_IGHJ.fasta

# #make blast database
# makeblastdb -in formatted_IGHV.fasta -dbtype nucl -parse_seqids -out IGHV
# makeblastdb -in formatted_IGHD.fasta -dbtype nucl -parse_seqids -out IGHD
# makeblastdb -in formatted_IGHJ.fasta -dbtype nucl -parse_seqids -out IGHJ

#######################################################################
#  filter IGHV transcripts and identify highly expressed transcripts  #
#######################################################################

printf "Filtering IGHV transcripts using blastn...\n\n"

#filter IGHV sequnces by blastn 
f=*Trinity.fasta
#get transcript ids having match with ig-genes by blastn
blastn -db ../data/IGHV -query $f -out $input_bam_file_name"_blastn" -num_threads $num_threads -outfmt 6;
#remove duplicated ids
less $input_bam_file_name"_blastn" | awk '!seen[$1]++' | awk '{print $1}' > $input_bam_file_name"_blastn_HG_IG_remdup.ID"
#filter transcript fasta file by transcript id
less $f | seqkit grep -f $input_bam_file_name"_blastn_HG_IG_remdup.ID" > $input_bam_file_name"_blastn_HG_IG_remdup.ID.fa"

printf "\n\n running salmon to estimate transcript abundace...\n\n"

#create sailfish index and quantify transcript abundances
fa=$input_bam_file_name"_blastn_HG_IG_remdup.ID.fa"
salmon index -t $fa -i salmon_index -k 25 -p $num_threads
salmon quant -i salmon_index -l "OSR" -1 $input_bam_file_name.slice.R1.fastq -2 $input_bam_file_name.slice.R2.fastq -o salmon_index

#sort transcript ids by tpm values
less ./salmon_index/quant.sf | sort -nr -k4 | grep -v "Name" | sed  '1i #Name    Length  EffectiveLength TPM     NumReads' >$input_bam_file_name.ig-transcripts.sortedbyTPM.txt
less $input_bam_file_name.ig-transcripts.sortedbyTPM.txt | awk '{print $1}' | grep -v "Name" >list

#sort fasta file by transcript ids
while read line; do
less $f | seqkit grep -p $line;
done <list >>$input_bam_file_name.ig-transcripts.sortedbyTPM.fasta

# igblastn for SHM and clonotypes
igblastn -germline_db_V ../data/IGHV -num_alignments_V 3 -germline_db_J ../data/IGHJ -num_alignments_J 3 -germline_db_D ../data/IGHD -num_alignments_D 3 -organism human -query $input_bam_file_name.ig-transcripts.sortedbyTPM.fasta -show_translation -auxiliary_data ../data/human_gl.aux >$input_bam_file_name.IgBLAST_out.txt

#copy all files to the outdir
cp $input_bam_file_name*txt $output_directory
cp $input_bam_file_name*fasta $output_directory

#remove temporary directory
rm -r ../CRIS_out.$input_bam_file_name

printf "#### finished job...\n\n"
