#!/usr/bin/bash
dir_genomes=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/05_GenomeAnnot/01_RepeatAnalysis/gordius_repeat
dir_rna=$HOME/data1/sjzhang/04_Nematomorpha/02_Transcriptome/01_Pre-assemble/01_fastp
dir_out=$HOME/data1/sjzhang/01_Nematodes/02_Rhabdias/01_Genome/05_GenomeAnnot/04_BrakerAnnot
if [[ ! -d ${dir_out} ]];then mkdir -p ${dir_out};fi
NP=48
NT=2

# tmp_fifofile="/tmp/$$.fifo"
# trap "exec 9>&-;exec 9<&-;exit 0" 2
# mkfifo ${tmp_fifofile}
# exec 9<>${tmp_fifofile}
# rm ${tmp_fifofile}
# for ((i=1;i<=${NT};i++));do echo >&9;done

for genomes in ${dir_genomes}/gordius.asm.FINAL.fasta.masked;do
	# read -u9
	# {

	out_name=$(basename ${genomes} .asm.FINAL.fasta.masked)
	out_dir=${dir_out}/${out_name}_braker2
	if [[ ! -d ${out_dir} ]];then mkdir -p ${out_dir};fi
	if [[ ! -d ${out_dir}/BUSCO_prot ]];then mkdir -p ${out_dir}/BUSCO_prot;fi
	cd ${out_dir} &&
	cp -r /usr/local/softwares/Augustus/config/ ${out_dir} &&
	echo -e "\n=============================\nBraker2 training for ${out_name}.\n$(date)\n=============================\n"
	
	braker.pl \
			  --species=Nematomorpha \
			  --genome=${genomes} \
			  --rnaseq_sets_ids=${out_name}-freeliving_clean,${out_name}-imature_clean \
			  --rnaseq_sets_dir=${dir_rna} \
			  --softmasking \
			  --AUGUSTUS_CONFIG_PATH=${out_dir}/config \
			  --gff3 &&
	echo -e "\n=============================\nbraker2 training for ${out_name} has been complted !\n$(date)\n=============================\n" ||
	echo -e "\n=============================\nCheck the setup .....\n=============================\n"
	echo -e "\n=============================\nStarting busco evaluation of braker2 predicted protein for ${out_name}.... !\n$(date)\n=============================\n"
	
	source /usr/local/anaconda3/bin/activate buscoenv
	busco_db=$HOME/data0/dbs/BUSCO_db/busco_downloads/lineages/metazoa_odb10
	cd ${busco_db} &&
	if [[ ! -f lengths_cutoff.bak ]];then cp lengths_cutoff lengths_cutoff.bak;fi
	if [[ ! -f scores_cutoff.bak ]];then cp scores_cutoff scores_cutoff.bak;fi
	lengths_cutoff_eff=2
	# scores_cutoff_eff=0.8
	# awk -v x=$scores_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; scr = $2; nscr = scr * x; print gen, nscr }' scores_cutoff.bak > scores_cutoff &&
	awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; B = $2; delta = $3; len =$4; ndelta = delta * x; print gen, B, ndelta, len }' lengths_cutoff.bak > lengths_cutoff 
	# awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; B = $2; delta = $3; len = $4; nlen = len * x; print gen, B, delta, nlen }' lengths_cutoff.bak > lengths_cutoff #if you use this line, the lengths_cutoff_eff should be in [0,1], recommended as 0.8.
	cd ${out_dir} &&
	busco -m prot \
		  -i ${out_dir}/braker/braker.aa \
		  --out_path ${out_dir}/BUSCO_prot \
		  -o ${out_name}_busco \
		  -c ${NP} \
		  --offline \
		  -l ${busco_db} 
	cp -f ${busco_db}/lengths_cutoff.bak ${busco_db}/lengths_cutoff &&
	cp -f ${busco_db}/scores_cutoff.bak ${busco_db}/scores_cutoff &&
	echo -e "\n=============================\nBusco evaluation of braker2 predicted protein for ${out_name} is completed !\n$(date)\n=============================\n"
	conda deactivate

	# echo >&9
	# } &
done
wait 

source /usr/local/anaconda3/bin/activate buscoenv
file_dir=${dir_out}/${out_name}_braker2/BUSCO_prot
cp_dir=${file_dir}/plots
if [[ ! -d ${cp_dir} ]];then mkdir -p ${cp_dir};fi
cd ${file_dir}
echo -e "========================================\nCopying the busco evaluation files ....\n========================================\n\n"

for dircts in ${file_dir}/*_busco;do
	if [[ -d ${dircts} ]];then
		species=$(basename ${dircts} _busco)
		echo -e "========================================\nCopying the busco evaluation files ${species}....\n========================================\n\n"
		echo -e "Copying ${species} file !!"
		cp -f ${dircts}/*.txt ${cp_dir}/
	else
		echo "${species} file is not a species !!"
	fi
done
echo -e "========================================\nCopying the busco evaluation files completed !!!\n========================================\n\n"
echo -e "========================================\nBusco ploting ....\n========================================\n\n"
cd ${cp_dir}

generate_plot.py -wd ${cp_dir}
echo -e "========================================\nBusco ploting completed !!\n========================================\n\n"

echo -e "\n=============================\nAll analysis have been completed !!\n$(date)\n=============================\n"





