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
hisat2 --dta \
       -p ${NP} \
       -x ${indexs}/${index_name}.fa \
       -1 ${out_name}_1.fq.gz \
       -2 ${out_name}_2.fq.gz \
       -S /dev/stdout 2>${out_name}_${index_name}.log | samtools view -bS - | samtools sort - > ${out_name}_${index_name}.sort.bam

# Merge the alignments and extract transcripts from aligments:
samtools merge -@ ${NP} -o merged.bam 1.bam 2.bam 3.ban ...
stringtie -p ${NP} -o stringtie.gtf merged.bam &&
gtf_genome_to_cdna_fasta.pl stringtie.gtf genome.fa.masked > transcripts.fasta &&
gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3

# Predict coding genes/peptieds using TransDecoder:
TransDecoder.LongOrfs -t transcripts.fasta -O ${out2dir} -m 100  --genetic_code universal

# Homology search
blastp -query *transdecoder_dir/longest_orfs.pep -db ${swissprot} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads ${NP} -out blastp.outfmt6 &
hmmscan --cpu ${NP} --domtblout pfam.domtblout ${Pfam} *transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --output_dir ${out2dir} ##this will generate "*transcripts.fasta.transdecoder.gff3"
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts_alignments.gff3 # merge the transcripts and protein coding gene annotations

```

