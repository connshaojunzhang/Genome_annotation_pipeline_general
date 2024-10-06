#!/usr/bin/bash
dir_genome=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/05_GenomeAnnot/01_RepeatAnalysis/gordius_autohic_post_repeat/gordius_autohic_post.fasta.masked
dir_protein=$HOME/data0/dbs/EnTAP_db/03_uniprot
dir_out=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/05_GenomeAnnot/05_HomologPredict
if [[ ! -d ${dir_out} ]];then mkdir -p ${dir_out};fi
NP=24
NT=2

# tmp_fifofile="/tmp/$$.fifo"
# trap "exec 9>&-;exec 9<&-;exit 0" 2
# mkfifo ${tmp_fifofile}
# exec 9<>${tmp_fifofile}
# rm ${tmp_fifofile}
# for ((i=1;i<=${NT};i++));do echo >&9;done

ln -s ${dir_genome} gordius_masked.fasta
for genome in ${dir_genome}/*.fasta;do


	out_name=$(basename ${genome} _masked.fasta)
	out_dir=${dir_out}/${out_name}_miniprot
	if [[ ! -d ${out_dir} ]];then mkdir -p ${out_dir};fi
	cd ${out_dir} &&
	echo -e "\n=============================\nBuilding miniprot index for genome of ${out_name}.\n$(date)\n=============================\n"
	
	miniprot -t${NP} -d genome.mpi ${genome} &&

	echo -e "\n=============================\nStarting Homology-based prediction for ${out_name} genome.... !\n$(date)\n=============================\n"

	miniprot -It${NP} --gff genome.mpi ${dir_protein}/uniprot_sprot.fasta > ${out_dir}/${out_name}_miniprot.gff3

done
	
echo -e "\n=============================\nAll analysis have been completed !!\n$(date)\n=============================\n"
