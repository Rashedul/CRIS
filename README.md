### CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq 

CRIS reconstructs entire *IGHV* gene to identify somatic hypermutation status in chronic lymphocytic leukemia using RNA-seq. CRIS has been validated against PCR-Sanger based clinical data for somatic hypermutations in CLL.

* Operating System: Linux

### Download the repository

```
git clone https://github.com/Rashedul/CRIS
cd CRIS/ # run CRIS from this directory 
```

### Dependencies

```
  - picard (v2.20.3)
  - trinity (v2.1.1)
  - blast (v2.9.0)
  - seqkit (v0.12.0)
  - sambamba (v0.7.0)
  - salmon (v0.8.1)
  - igblast (v1.14.0)
```

Executables must be accessible from user's PATH. CRIS has been tested on the software versions indicated.

### Installation of dependencies using conda

```
conda env create --file environment.yml --force
conda activate cris_env
```

### Usage

```
bash CRIS.sh -h

Construction of IHGV transcripts from RNA-seq.

Usage: CRIS.sh -inbam <input_bam_file> -threads <num_threads> -memory <max_memory_assembly>
                        <input_bam_file>: bam file must be aligned to hg38 genome build, coordinate-sorted and indexed
                        <num_threads>: number of threads; default 4
                        <max_memory_assembly>: maximum memory in G (gigabyte) allowed for assembly; default 4G
```

### Test run

**NOTE:** bam file must be aligned to the hg38 genome build, coordinate-sorted and indexed. This pipeline is tested for BWA and STAR alignment of paired-end RNA-seq reads. 
 
IGHV status for SRR1814049 (US-1422278) using Sanger sequencing:

* IGHV gene: IGHV3-74
* Percent identity: 94.6

```
bash CRIS.sh -inbam SRR1814049_test.bam
or 
bash CRIS.sh -inbam /path/SRR1814049_test.bam
```

### Output files

```
`SRR1814049_test.bam.ig-transcripts.sortedbyTPM.fasta` contains Ig-transcript fasta sequences ordered by expression (TPM) values. 

`SRR1814049_test.bam.ig-transcripts.sortedbyTPM.txt` contains Ig-transcript sequence IDs and expression values.

`SRR1814049_test.bam.IgBLAST_out.txt` contains the percent identity and alignment between Ig-transcript and top germline V gene hit.
```

### License 

This project is licensed under the [MIT license](https://github.com/Rashedul/CRIS/blob/main/LICENSE).

### Contact 

Rashedul Islam (rashed1 (at) student.ubc.ca)

