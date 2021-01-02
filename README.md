### CRIS: Complete Reconstruction of Immunoglobulin V-D-J Sequences from RNA-seq 

Construction of IHGV transcripts from RNA-seq

### Installation

```
git clone https://github.com/Rashedul/CRIS
cd CRIS

conda env create --file environment.yml --force
source activate cris_env

```

### Test run

Sanger output for SRR1814049 (US-1422278): IGHV gene: V3-74; percent identity: 94.6

```
bash CRIS.sh -h

bash CRIS.sh -inbam ./test_data/SRR1814049_subset.bam

```


### Output 

```
`file.bam.ig-transcripts.sortedbyTPM.fasta` contains Ig fasta sequences ordered by expression (TPM) values.

`file.bam.ig-transcripts.sortedbyTPM.fasta` contains Ig sequences IDs and expression values.

```

