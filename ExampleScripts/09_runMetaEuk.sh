#!/usr/bin/bash
dir_genomes=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/05_GenomeAnnot/01_RepeatAnalysis/gordius_repeat
dir_out=$HOME/data1/sjzhang/04_Nematomorpha/01_Genome/05_GenomeAnnot/03_Metaeuk_denovo
if [[ ! -d ${dir_out} ]];then mkdir -p ${dir_out};fi

metaeuk_db=$HOME/data0/dbs/metaeukdb

source /usr/local/anaconda3/bin/activate annotations

NP=48
NT=2

tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo ${tmp_fifofile}
exec 9<>${tmp_fifofile}
rm ${tmp_fifofile}
for ((i=1;i<=${NT};i++));do echo >&9;done
cd ${dir_out} &&
ln -s ${dir_genomes}/gordius.asm.FINAL.fasta.masked gordius.masked.fasta

for dbs in ${metaeuk_db}/*UniRef*;do
	read -u9
	{

	dbs_name=$(basename ${dbs})
	out_dir=${dir_out}/gordius_${dbs_name#*_}
	if [[ ! -d ${out_dir} ]];then mkdir -p ${out_dir};fi

	echo -e "\n=============================\nMetaEuk prediction for Gordius is started.\n$(date)\n=============================\n"
	metaeuk easy-predict gordius.masked.fasta ${dbs}/*.metaeukdb ${out_dir}/gordius_${dbs_name#*_} ${out_dir}/tmp \
			--min-length 15 \
			-e 100 \
			--metaeuk-eval 0.0001 \
			--metaeuk-tcov 0.5 \
			--exhaustive-search 1 \
			--num-iterations 1 \
			--mask-lower-case 0 \
			--threads ${NP} &&
	echo -e "\n=============================\nMetaEuk prediction for Gordius has been complted !\n$(date)\n=============================\n" ||
	echo -e "\n=============================\nCheck the setup .....\n=============================\n"

	echo >&9
	} &
done
wait 
echo -e "\n=============================\nAll analysis have been completed !!\n$(date)\n=============================\n"
conda deactivate

busco_db=$HOME/data0/dbs/BUSCO_db/busco_downloads/lineages/metazoa_odb10
##activate the envs of busco and quanst
source /usr/local/anaconda3/bin/activate buscoenv

NP=20
NM=100
# NT=$[$(ls -lh ${genome_dir}/*_genome/*.gz | wc -l)]
# NT=3

echo -e echo -e "=========================================\nChange the BUSCO config file !\n=========================================\n\n"
cd ${busco_db} 
if [[ ! -f lengths_cutoff.bak ]];then cp lengths_cutoff lengths_cutoff.bak;fi
if [[ ! -f scores_cutoff.bak ]];then cp scores_cutoff scores_cutoff.bak;fi
lengths_cutoff_eff=2
# scores_cutoff_eff=0.8
# awk -v x=$scores_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; scr = $2; nscr = scr * x; print gen, nscr }' scores_cutoff.bak > scores_cutoff &&
awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; B = $2; delta = $3; len = $4; ndelta = delta * x; print gen, B, ndelta, len }' lengths_cutoff.bak > lengths_cutoff 
# awk -v x=$lengths_cutoff_eff 'BEGIN{OFS=FS="\t"}{gen = $1; B = $2; delta = $3; len = $4; nlen = len * x; print gen, B, delta, nlen }' lengths_cutoff.bak > lengths_cutoff #if you use this line, the lengths_cutoff_eff should be in [0,1], recommended as 0.8.

tmp_fifofile="/tmp/$$.fifo"
trap "exec 9>&-;exec 9<&-;exit 0" 2
mkfifo ${tmp_fifofile}
exec 9<>${tmp_fifofile}
rm ${tmp_fifofile}
for ((i=1;i<=${NT};i++));do echo >&9;done

##geting the statistics for the assemblies

# if [[ ! -d ${dir_out}/01_quast_Stats ]];then mkdir -p ${dir_out}/01_quast_Stats;fi
# echo -e "Assembley\tcontigs\tlength\tGC\tN50\tL50" > ${dir_out}/01_quast_Stats/assemble.statistics.txt
cd ${dir_out} &&
for genomes in ${dir_out}/gordius*;do
	read -u9
	{

	out_name=$(basename ${genomes})
	# out2dir=${dir_out}/01_quast_Stats/${out_name}_quasta
	# if [[ ! -d ${out2dir} ]];then mkdir -p ${out2dir};fi
	# # echo -e "=========================================\nStatistics of N50 after spades assembly of ${out_name} started:\n$(date)\n=========================================\n\n"
	# quast -o ${out2dir} -t ${NP} ${genomes}  && 
	# cd ${out2dir} && sed -n "/${out_name}/p" transposed_report.txt | awk '{print $1"\t"$14"\t"$16"\t"$17"\t"$18"\t"$21}' >> ${dir_out}/assemble.statistics.txt &&
	# echo -e "=========================================\nStatistics of N50 is saved in ${dir_out}/assemble.statistics.txt\n$(date)\n=========================================\n\n"
	echo -e "=========================================\nBUSCO evaluation for ${ou_tname} is started at:\n$(date)\n=========================================\n\n"
	# gunzip -c ${genomes} > ${genomes%.gz} &&
	out_put=${dir_out}/02_BUSCO_prot
	if [[ ! -d ${out_put} ]];then mkdir -p ${out_put};fi
	busco -m prot \
		  -i ${genomes}/gordius_Uni*0.fas \
		  --out_path ${out_put} \
		  -o ${out_name}_busco \
		  -c ${NP} \
		  --offline \
		  -l ${busco_db} &&
	echo -e "=========================================\nBUSCO evaluation for ${out_name} is completed !\n$(date)\n=========================================\n\n" || echo -e "=========================================\n BUSCO evaluation for ${outname} has some problems !! \n=========================================\n\n"

	echo >&9
	} &
done
wait

cp -f ${busco_db}/lengths_cutoff.bak ${busco_db}/lengths_cutoff &&
cp -f ${busco_db}/scores_cutoff.bak ${busco_db}/scores_cutoff &&
echo -e "=========================================\n Statistics of N50 and BUSCO after assembling are all completed, ploting BUSCO now !!!!\n$(date)\n=========================================\n\n"

# ploting the BUSCO evaluation
file_dir=${dir_out}/02_BUSCO_prot
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
conda deactivate




