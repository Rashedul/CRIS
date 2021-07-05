### CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq 

CRIS reconstructs entire *IGHV* gene to identify somatic hypermutation status in chronic lymphocytic leukemia using RNA-seq. CRIS has been validated against PCR-Sanger based clinical data for somatic hypermutations in CLL.

### CRIS in Docker container

* Minimum 16GB memory and 4 CPUs are required to be used by the container

Build and run docker image

```
# download
git clone https://github.com/Rashedul/CRIS
cd CRIS/ # run CRIS from this directory 

# build docker image
sudo docker build --tag cris:v1 . 

# run CRIS using test bam file. bam file must be aligned to the hg38 genome build and coordinate-sorted.
sudo docker run -v $PWD/SRR1814049_test.bam:/cris/SRR1814049_test.bam -t cris:v1 bash CRIS_docker.sh SRR1814049_test.bam 4 16G

```

### CRIS in bash 

* Operating System: Linux
* Dependencies: ```picard (v2.20.3), trinity (v2.1.1), blast (v2.9.0), seqkit (v0.12.0), sambamba (v0.7.0), salmon (v0.8.1), igblast (v1.14.0), jellyfish (v2.2.10)```. Executables must be accessible from user's PATH. 
* Install dependencies using conda environment

```
# download
git clone https://github.com/Rashedul/CRIS
cd CRIS/ # run CRIS from this directory 

# installing dependencies using conda
conda create --name cris_env --file environment.txt
conda activate cris_env 

# usage
bash CRIS.sh -h

CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq.

Usage: CRIS.sh -inbam <input_bam_file> -outdir <output_directory> -threads <num_threads> -memory <max_memory_assembly>
                        <input_bam_file>: (required) bam file must be aligned to hg38 genome build, coordinate-sorted and indexed
                        <output_directory>: (optional) full path of output directory or output files will be written in current directory
                        <num_threads>: (optional) number of threads; default 4
                        <max_memory_assembly>: (optional) maximum memory in G (gigabyte) allowed for assembly; default 16G

# test run from CRIS directory. bam file must be aligned to the hg38 genome build, coordinate-sorted and indexed.
bash CRIS.sh -inbam SRR1814049_test.bam 
or,
bash CRIS.sh -inbam /fullPath/SRR1814049_test.bam -outdir /fullPath/
```

### Output 

```
# expected mutational status for test run 

- IGHV gene: IGHV3-74
- Percent identity: 94.6

# output files

`SRR1814049_test.bam.IgBLAST_out.txt` contains the percent identity and alignment between Ig-transcript and top germline V gene hits.

`SRR1814049_test.bam.ig-transcripts.sortedbyTPM.fasta` contains Ig-transcript fasta sequences ordered by expression (TPM) values. 

`SRR1814049_test.bam.ig-transcripts.sortedbyTPM.txt` contains Ig-transcript sequence IDs and expression values.
```

### License 

This project is licensed under the [MIT license](https://github.com/Rashedul/CRIS/blob/main/LICENSE).

### Contact 

Rashedul Islam (rashed1 (at) student.ubc.ca)

