#!/bin/bash 
if [ "$1" == "-h" ] ; then
        echo -e "Construction of IHGV transcripts from RNA-seq.

Usage: `basename $0` -inbam <input_bam_file> -threads <num_threads> -memory <max_memory_assembly>
                        <input_bam_file>: bam file must be aligned to hg38 genome build, coordinate-sorted and indexed
                        <num_threads>: number of threads; default 4
                        <max_memory_assembly>: maximum memory in G (gigabyte) allowed for assembly; default 4G"
exit 0
fi

num_threads=4
max_memory_assembly=4G

while [ $# -gt 0 ]
do
        case "$1" in
        -inbam) input_bam_file="$2"; shift;;
        -threads) num_threads="$2"; shift;;
        -memory) max_memory_assembly="$2"; shift;;
        esac
        shift
done 

###############################################################
# Slice bam; PE bam to fastq; de novo assembly of transcripts #
###############################################################

mkdir CRIS.$input_bam_file
cd CRIS.$input_bam_file

#hg38 ig coordinates
echo "chr14   105600001       106880000
chr14_KI270726v1_random 0       43739
chr15   21710000        22190000
chr16   31950001        33970000
chr16_KI270728v1_random 0       1872759" >Ig_regions.bed

#slice bam by regions of interest. For slice, input file must be coordinate-sorted and indexed
sambamba slice ../$input_bam_file -L Ig_regions.bed -o $input_bam_file.slice.bam
picard SamToFastq I=$input_bam_file.slice.bam F=$input_bam_file.slice.R1.fastq F2=$input_bam_file.slice.R2.fastq

printf "Trinity is constructing de novo transcripts...\n\n"

#trinity
Trinity --seqType fq --max_memory $max_memory_assembly --left $input_bam_file.slice.R1.fastq --right $input_bam_file.slice.R2.fastq --CPU $num_threads --trimmomatic --full_cleanup;

#######################################################################
#  filter IGHV transcripts and identify highly expressed transcripts  #
#######################################################################

printf "Filtering IGHV transcripts using blastn...\n\n"

# download germline IGHV sequences
wget http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta
curl -O http://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta

#replace periods with Ns
less IGHV.fasta | sed -e '/^[^>]/s/[^ATGCatgc]/N/g' >IGHV_N.fasta

#make balst database
makeblastdb -in IGHV_N.fasta -dbtype nucl -parse_seqids -out IGHV

#filter IGHV sequnces by blastn 
f=*Trinity.fasta
#get transcript ids having match with ig-genes by blastn
blastn -db IGHV -query $f -out $input_bam_file"_blastn" -num_threads $num_threads -outfmt 6;
#remove duplicated ids
less $input_bam_file"_blastn" | awk '!seen[$1]++' | awk '{print $1}' > $input_bam_file"_blastn_HG_IG_remdup.ID"
#filter transcript fasta file by transcript id
less $f | seqkit grep -f $input_bam_file"_blastn_HG_IG_remdup.ID" > $input_bam_file"_blastn_HG_IG_remdup.ID.fa"

printf "running sailfish to estimate transcript abundace...\n\n"

#create sailfish index 
fa=$input_bam_file"_blastn_HG_IG_remdup.ID.fa"
sailfish index -t $fa -o sailfish_index -k 25 -p $num_threads

sailfish quant -i sailfish_index -l "OSR"  -1 $input_bam_file.slice.R1.fastq -2 $input_bam_file.slice.R2.fastq -o sailfish_index -p $num_threads

#sort transcript ids by tpm values
less ./sailfish_index/quant.sf | sort -nr -k4 | grep -v "Name" | sed  '1i #Name    Length  EffectiveLength TPM     NumReads' >$input_bam_file.ig-transcripts.sortedbyTPM.txt

less $input_bam_file.ig-transcripts.sortedbyTPM.txt | awk '{print $1}' | grep -v "Name" >list

#sort fasta file by transcript ids
while read line; do
less $f | seqkit grep -p $line;
done <list >>$input_bam_file.ig-transcripts.sortedbyTPM.fasta

#remove temporary files
rm -r sailfish_index
rm IGHV* *slice* *_blastn* *.Trinity.fasta *bed list

printf "#### finished job...\n\n"
