#!/bin/bash 
# cat "${1}"

# cat "${1}" | awk '$1>2{print}' >new.txt

# echo "print new file"
# cat new.txt

#hg38 ig coordinates
echo "chr14   105600001       106880000
chr14_KI270726v1_random 0       43739
chr15   21710000        22190000
chr16   31950001        33970000
chr16_KI270728v1_random 0       1872759" >Ig_regions.bed

echo "print ig regions"
cat Ig_regions.bed

#slice bam by regions of interest. For slice, input file must be coordinate-sorted and indexed
#sambamba -h

input_bam_file_name="${1}"
ls -l $input_bam_file_name
echo $input_bam_file_name.slice.bam

sambamba index $input_bam_file_name
sambamba slice $input_bam_file_name -L Ig_regions.bed -o $input_bam_file_name.slice.bam
sambamba flagstat $input_bam_file_name.slice.bam

picard SamToFastq I=$input_bam_file_name.slice.bam F=$input_bam_file_name.slice.R1.fastq F2=$input_bam_file_name.slice.R2.fastq

ls -l *

#cat /etc/*ease