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
RepeatMasker -pa 20 -a -xsmall -s -gff -no_is -lib db_name-families.fa \
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
samtools merge -@ ${NP} -o merged.bam 1.bam 2.bam 3.bam ...
stringtie -p ${NP} -o stringtie.gtf merged.bam &&
gtf_genome_to_cdna_fasta.pl stringtie.gtf genome.fa.masked > transcripts.fasta &&
gtf_to_alignment_gff3.pl stringtie.gtf > transcripts.gff3

# Predict coding genes/peptieds using TransDecoder:
TransDecoder.LongOrfs -t transcripts.fasta -O ${out2dir} -m 100  --genetic_code universal

# Homology search
blastp -query *transdecoder_dir/longest_orfs.pep -db ${swissprot} -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads ${NP} -out blastp.outfmt6 &
hmmscan --cpu ${NP} --domtblout pfam.domtblout ${Pfam} *transdecoder_dir/longest_orfs.pep
TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits pfam.domtblout --retain_blastp_hits blastp.outfmt6 --output_dir ${out2dir} ##this will generate "*transcripts.fasta.transdecoder.gff3"
cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts_transdecoder.gff3 # merge the transcripts and protein coding gene annotations
# transcripts_transdecoder.gff3 will be generated.
```
#### 2. Homology-based gene prediciton (protein alignments with splicing)
* genome.fa.masked
* protein database (e.g., homology proteins OrthoDB Metazoa.fa, or a more comprehesive protein database Uniprot/Swiss-prot)
```
# Build index for the genomes
miniprot -t${NP} -d genome.mpi genome.fa.masked
# Search
miniprot -It${NP} --gff genome.mpi uniprot_sprot.fasta > miniprot.gff3

# Reform the miniprot.gff3 to be compatible with downstream EVM pipline (this script ${after_miniprot_GFF_2_EVM_alignment_GFF3} can be found in ExampleScripts in this repository)
python ${after_miniprot_GFF_2_EVM_alignment_GFF3} miniprot.gff3 > protein_alignments.gff3
```
#### 3. Ab initio annotation
* genome.fa.masked
* proteins from OrthoDB, e.g., Metazoa.fa
* alignments of RNA-seqs, e.g., merged.bam used in "1. Gene prediction using RNA-seq"
```
# BRAKER2 pipline is used (including the GeneMark-ETP and augustus training)
braker.pl \
     --species=Your_species \ 
     --genome=genome.fa.masked \
     --prot_seq=${dir_orthoprots}/Metazoa.fa \
     --bam=merged.bam \
     --AUGUSTUS_CONFIG_PATH=${out_dir}/config \   # using the augustus config file.
     --gff3 \
     --threads ${NP}
# Rename the generated GFF3 file to be compatible with EVM pipline (${gff_rename} can be found in ExampleScripts in this repository)
python ${gff_rename} braker.gff3 Your_species > geneprediction_braker.gff3
```
#### 4. Gene structure annotation based on _de novo_ assembled transcripts
* genome.fa.masked
* RNA-seqs (e.g., for different develop stage or different timepoint)
##### 1. Transcriptome assemble using ORP or trinity
```
Trinity --seqType fq \
        --max_memory 500G --CPU ${NP} \
        --left RNA1_1.fq.gz,RNA2_1.fq.gz \
        --right RNA1_2.fq.gz,RNA2_2.fq.gz \
        --output ${prefix}
# ${prefix}.Trinity.fasta will be generated. 
```
##### 2. Gene structure annotation using PASA pipeline
```
# config the pasa database
cp $PASAHOME/pasa_conf/pasa.alignAssembly.Template.txt ${your_work_directory}/alignAssembly.config
cd ${your_work_directory}
sed -i "s#<__DATABASE__>#${your_work_directory}/${prefix}_pasadb#g" alignAssembly.config # provide full path to pasa database.

# clean the assembled transcript
seqclean ${prefix}.Trinity.fasta

# run pasa pipeline
$PASAHOME/Launch_PASA_pipeline.pl \
			-c alignAssembly.config -C -R \
			-g genome.fa.masked \
			-t ${prefix}.Trinity.fasta.clean -T \
			-u ${prefix}.Trinity.fasta \
			--ALIGNERS blat,gmap,minimap2
# this will generate ${prefix}_pasadb.pasa_assemblies.gff3, simply merge this annotation with GFF3 generated by transdecoder (transcripts_transdecoder.gff3).
cat ${prefix}_pasadb.pasa_assemblies.gff3 transcripts_transdecoder.gff3 > transcripts_alignments.gff3
```
## Step3: Merge the annotation from RNA-seq, homology-based, ab initio strategies
* genome.fa.masked
* transcripts_alignments.gff3, protein_alignments.gff3, and geneprediction_braker.gff3
* weight file from EVM pipline
```
# generate weight.txt
echo -e "ABINITIO_PREDICTION\tEVM\t10" >> weights.txt
echo -e "PROTEIN\tminiprot_protAln\t5" >> weights.txt
echo -e "TRANSCRIPT\tassembler-${prefix}_pasadb\t10" >> weights.txt
echo -e "OTHER_PREDICTION\ttransdecoder\t10" >> weights.txt

EVidenceModeler \
   --sample_id ${prefix} \
   --weights weights.txt \
   --genome genome_masked.fasta \
   --gene_predictions gene_prediction.gff3 \
   --protein_alignments protein_alignments.gff3 \
   --transcript_alignments transcript_alignments.gff3 \
   --segmentSize 100000 \
   --overlapSize 10000 \
   --debug
# this will generate ${prefix}.EVM.gff3
```
## Step4: For ncRNA, including lncRNA annotation
