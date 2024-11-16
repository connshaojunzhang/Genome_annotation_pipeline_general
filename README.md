# Genome_annotation_pipeline_general (_de novo_ assembled genome )

## Data:
  * genome.fa (chromosome level assembly, e.g., after AutoHiC assisted assemble).
  * RNA-seqs (with different development stage or from different organs, same speices of the  genome.fa).
## Step1: Repeat elements annotation and masking.
  **1.RepeatModeler & RepeatMasker**
  * genome.fa  
#### Building database of repeat elements
```
BuildDatabase -name db_name genome.fa
```
#### Predict repeat elements using RepeatModeler
```
RepeatModeler -database db_name -threads 20 -LTRStruct
```
#### Genome masking using RepeatMasker
```
RepeatMasker -pa 20 -a -s -gff -no_is -lib db_name-families.fa \
             -dir outdir genome.fa
```
* genome.fa.masked (see genome.fa.tbl for repeat elements statistics)
#### For transposon elements annotation
```
EDTA.pl --genome genome.fa --species others --sensitive 1 --anno 1 --evaluate 1 --threads 24
```
## Step2: Gene prediction and annotation.
#### 1. Gene prediction using RNA-seq
* genome.fa.masked
* RNA-seqs

```
# Transcripts alignment to genome using hisat2
hisat2-build -p 24 genome.fa.masked genome_index


```

