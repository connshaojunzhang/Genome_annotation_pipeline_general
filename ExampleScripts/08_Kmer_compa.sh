#!/usr/bin/bash
dir_seq=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/01_Assemble3rd/01_Pre-assemble/02_deconta/Gordius0725_deconta
dir_genomes=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/04_Hic_assemble/06_3D-DNA_post
dir_out=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/05_GenomeAnnot/02_Evaluation/03_KmerComp
if [[ ! -d ${dir_out} ]];then mkdir -p ${dir_out};fi

NP=64
NT=2

tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo ${tmp_fifofile}
exec 9<>${tmp_fifofile}
rm ${tmp_fifofile}
for ((i=1;i<=${NT};i++));do echo >&9;done

source /usr/local/anaconda3/bin/activate kmerenv

for genomes in ${dir_genomes}/gordius.asm.FINAL.fasta;do
	# read -u9
	# {
	cd ${dir_out} 

	out_name=$(basename ${genomes} .asm.FINAL.fasta)
	for kmer in $(seq 23 4 43);do
		read -u9
		{

		out_put=${dir_out}/${out_name}.kcomp.${kmer}mer
		kat comp -t ${NP} -m ${kmer} -n -p pdf -h \
				 -o ${out_put} \
				 ${dir_seq}/Gordius0725_deconta_1.fq.gz \
				 ${dir_seq}/Gordius0725_deconta_2.fq.gz \
				 ${genomes} &&
		echo -e "\n=============================\nkat_comp for ${out_name}_${kmer}mer has been completed !!\n$(date)\n=============================\n" ||
		echo -e "\n=============================\nkat_comp for ${out_name}_${kmer}mer has some problems!\n$(date)\n=============================\n"

		echo >&9
		} &
	done
	wait
	for kmer in $(seq 23 4 43);do
		kat plot spectra-cn -o ${out_name}.kcomp.${kmer}mer-main -x 200 -y 2000000 ${out_name}.kcomp.${kmer}mer-main.mx &&
		echo -e "\n=============================\nkat_plot for ${out_name}_${kmer}mer has been completed !!\n$(date)\n=============================\n" ||
		echo -e "\n=============================\nkat_plot for ${out_name}_${kmer}mer has some problems!\n$(date)\n=============================\n"
	done
done

conda deactivate
echo -e "\n=============================\n===All analysis have been completed !!===\n=============================\n"




