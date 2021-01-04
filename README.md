### CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq 

Construction of IHGV transcripts from RNA-seq

### Installation

* Operating System: Linux

```
git clone https://github.com/Rashedul/CRIS
cd CRIS

conda env create --file environment.yml --force
source activate cris_env
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

**NOTE:** bam file must be aligned to the hg38 genome build, coordinate-sorted and indexed. This pipeline is tested on BWA alignment of RNA-seq reads to the hg38 genome build.
 
IGHV status for SRR1814049 (US-1422278) using Sanger sequencing:

* IGHV gene: IGHV3-74
* Percent identity: 94.6

```
bash CRIS.sh -inbam SRR1814049_test.bam
```

### Output 

Upload the Ig fasta file to [IgBlast](https://www.ncbi.nlm.nih.gov/igblast/) to get the percent identity between query and top germline V gene hit.

```
`SRR1814049_test.bam.ig-transcripts.sortedbyTPM.fasta` contains Ig fasta sequences ordered by expression (TPM) values. Upload this fasta file to IgBlast.

`SRR1814049_test.bam.ig-transcripts.sortedbyTPM.txt` contains Ig sequence IDs and expression values.
```

### License 

This project is licensed under the [MIT license](https://github.com/Rashedul/CRIS/blob/main/LICENSE).
