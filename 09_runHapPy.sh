#!/usr/bin/bash
dir_seq=$HOME/data0/Genomes/sjzhang/02_Nematomorpha/01_PacbioHiFI
dir_genomes=$HOME/data3/sjzhang/02_Nematomorpha/01_Genome/04_Hic_assemble/06_3D-DNA_post
dir_out=$HOME/data3/sjzhang/02_Nematomorpha/01_Genome/05_GenomeAnnot/02_Evaluation/04_HaploEstimate
if [[ ! -d ${dir_out} ]];then mkdir -p ${dir_out};fi

NP=64
NT=2

# tmp_fifofile="/tmp/$$.fifo"
# trap "exec 9>&-;exec 9<&-;exit 0" 2
# mkfifo ${tmp_fifofile}
# exec 9<>${tmp_fifofile}
# rm ${tmp_fifofile}
# for ((i=1;i<=${NT};i++));do echo >&9;done

source /usr/local/anaconda3/bin/activate kmerenv

for genomes in ${dir_genomes}/gordius.asm.FINAL.fasta;do
	# read -u9
	# {
	cd ${dir_out} 

	out_name=$(basename ${genomes} .asm.FINAL.fasta)
	out_put=${dir_out}/${out_name}_haploEstimate
	if [[ ! -d ${out_put} ]];then mkdir -p ${out_put};fi
	cd ${out_put} &&
	minimap2 -t ${NP} -ax map-pb ${genomes} ${dir_seq}/gordius.hifi_reads.fa --secondary=no | \
	samtools sort - ${out_name}.sort.bam &&
	samtools index ${out_name}.sort.bam &&
	Hap.py coverage -t ${NP} -d ./ ${out_name}.sort.bam &&
	Hap.py estimate -S 230M -p ${out_name}.Hap.estimate --plot ./${out_name}.sort.bam.hist &&
	echo -e "\n=============================\nHapPy for ${out_name}_${kmer}mer has been completed !!\n$(date)\n=============================\n" ||
	echo -e "\n=============================\nHapPy for ${out_name}_${kmer}mer has some problems!\n$(date)\n=============================\n"
python /usr/local/softwares/HapPy-0.1/happy/Hap.py estimate -ll 1 -ld 50 -lh 120 -s 230M -o ./happy_stats --plot ./happy_out/gordius.sort.bam.hist 
		# echo >&9
		# } &
done
	
conda deactivate
echo -e "\n=============================\n===All analysis have been completed !!===\n=============================\n"




