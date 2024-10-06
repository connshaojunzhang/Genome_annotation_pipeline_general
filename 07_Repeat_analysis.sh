#!/usr/bin/bash
dir_genomes=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/04_Hic_assemble/06_3D-DNA_post
dir_out=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/05_GenomeAnnot/01_RepeatAnalysis
if [[ ! -d ${dir_out} ]];then mkdir -p ${dir_out};fi
repeatmodeler=dfam/tetools:latest
docker_mount=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome
users=$(echo $USER)
# NT=2

# tmp_fifofile="/tmp/$$.fifo"
# trap "exec 9>&-;exec 9<&-;exit 0" 2
# mkfifo ${tmp_fifofile}
# exec 9<>${tmp_fifofile}
# rm ${tmp_fifofile}
# for ((i=1;i<=${NT};i++));do echo >&9;done

for genomes in ${dir_genomes}/gordius.asm.FINAL.fasta;do
	read -u9
	{

	out_name=$(basename ${genomes} .asm.FINAL.fasta)
	docker_home=/home/docker
	out_dir=${dir_out}/${out_name}_repeat
	if [[ ! -d ${out_dir} ]];then mkdir -p ${out_dir};fi
	echo -e "\n=============================\nStep 01: Building repeat database for ${out_name}.\n$(date)\n=============================\n"
	docker run --user=root --rm -v ${docker_mount}:${docker_home} ${repeatmodeler} \
			BuildDatabase -name ${docker_home}/05_GenomeAnnot/01_RepeatAnalysis/${out_name}_repeat/${out_name} ${docker_home}/04_Hic_assemble/06_3D-DNA_post/${out_name}.asm.FINAL.fasta &&
	echo -e "\n=============================\nStep 02: Starting RepeatModeler for ${out_name}.\n$(date)\n=============================\n"
	sudo chown -R ${users}:${users} ${out_dir} &&
	docker run --user=root --rm -v ${docker_mount}:${docker_home} -w ${docker_home}/05_GenomeAnnot/01_RepeatAnalysis/${out_name}_repeat ${repeatmodeler} \
			RepeatModeler -database ${docker_home}/05_GenomeAnnot/01_RepeatAnalysis/${out_name}_repeat/${out_name} \
						  -threads 20 -LTRStruct &&
	sudo chown -R ${users}:${users} ${out_dir} &&
	echo -e "\n=============================\nRepeatModeler for ${out_name} completed.\n$(date)\nStarting Step03: RepeatMasker for ${out_name} now.\n================================\n"
	docker run --user=root --rm -v ${docker_mount}:${docker_home} -w ${docker_home}/05_GenomeAnnot/01_RepeatAnalysis/${out_name}_repeat ${repeatmodeler} \
			RepeatMasker -pa 20 -a -s -gff -no_is -lib ${docker_home}/05_GenomeAnnot/01_RepeatAnalysis/${out_name}_repeat/${out_name}-families.fa -dir ${docker_home}/05_GenomeAnnot/01_RepeatAnalysis/${out_name}_repeat ${docker_home}/04_Hic_assemble/06_3D-DNA_post/${out_name}.asm.FINAL.fasta &&
	sudo chown -R ${users}:${users} ${out_dir} &&
	echo -e "\n=============================\nRepeatModeler for ${out_name} has been completed !!\n$(date)\n=============================\n"

	echo >&9
	} &
done
wait 
echo -e "\n=============================\n===All analysis have been completed !!===\n=============================\n"




