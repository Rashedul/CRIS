### CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq 

CRIS reconstructs the immunoglobulin heavy chain variable region *IGHV* gene, enumerates single nucleotide variants and predicts hypermutation status from RNA-seq datasets. Both ribodepleted and polyA selected RNA-seq datasets are appropriate with a minimum of 25M sequence reads per sample. CRIS has been validated against clinical PCR-Sanger based hypermutation classification in the context of Chronic Lymphocytic Leukemia.

### Requirements

* A minimum of 16GB RAM and 4 threads are required to run CRIS 
* bam file must be aligned to hg38 (GRCh38) genome build using BWA, coordinate-sorted and indexed

### CRIS in [Docker](https://docs.docker.com/) container 

Build and run docker image. 

```
# download
git clone https://github.com/Rashedul/CRIS
cd CRIS/ # run CRIS from this directory 

# build docker image
sudo docker build --tag cris:v1 . 

# to confirm install please run CRIS using the supplied test bam file (SRR1814049_test.bam) that has been aligned to hg38 (GRCh38) build and coordinate-sorted using SAMBAMBA sort. As mentioned 4 threads and 16G RAM is allocated by CRIS by default.
sudo docker run --name cris_analysis -v /fullPath/SRR1814049_test.bam:/cris/SRR1814049_test.bam -t cris:v1 bash CRIS_docker.sh SRR1814049_test.bam threads=4 memory=16G 

# copy output files to local path
sudo docker cp cris_analysis:cris/ .

# export docker container to local path, output is stored within cris directory.
sudo docker export cris_analysis > cris_analysis_container.tar

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

* [Example output files are here.](https://github.com/Rashedul/CRIS/tree/main/data/output_examples)


```
# expected mutational status for test run 

- IGHV gene: IGHV3-74
- Percent identity: 94.6

# output files

- SRR1814049_test.bam.IgBLAST_out.txt
Contains the percent identity and alignment between Ig-transcript and top germline IGHV gene hits.

- SRR1814049_test.bam.ig-transcripts.sortedbyTPM.txt 
Contains Ig-transcript sequence IDs and sorted by expression (TPM) values.

- SRR1814049_test.bam.ig-transcripts.sortedbyTPM.fasta 
Contains Ig-transcript fasta sequences sorted by expression (TPM) values. 

```

### License 

This project is licensed under the [MIT license](https://github.com/Rashedul/CRIS/blob/main/LICENSE).

### Contact 

Rashedul Islam (rashedul.islam@alumni.ubc.ca)


### Publication 

Rashedul Islam, Misha Bilenky, Andrew P. Weng, Joseph M. Connors, Martin Hirst. [CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq data](https://academic.oup.com/bioinformaticsadvances/advance-article/doi/10.1093/bioadv/vbab021/6367791?login=true). *Bioinformatics Advances*, vbab021, September 2021.