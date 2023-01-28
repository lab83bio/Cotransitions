#! /usr/bin/env bash
#cotr pipeline

##################
# DEFINE VARIABLES
##################

level="Eukaryota" #"Bacteria", "Archaea" (faster), "Mammalia", etc
taxonomy_database='ncbi' #ncbi|oma
tree="raxml" #raxml|ncbi|random
ladder=("RL" "LL" "NL") #tree orientation (RL=right-ladderized)
ncores=30

#################################################
# MAKE OUTDIR (example: "Eukaryota_20-09-2022_0")
#################################################

cwd=`realpath .`
date=`date '+%d-%m-%Y'`
n=1
while ! mkdir ${level}_HOGs_${date}_${n}
do
    n=$((n+1))
done
cd ${level}_HOGs_${date}_${n}

#################
# START PIPLELINE
#################

curl https://v101.orthodb.org/download/odb10v1_levels.tab.gz -o HOGs_${level}_${taxonomy_database}_levels.tab.gz 2>>log.txt
${cwd}/Utilities/OMA_translator.py -t ${level} -d ${taxonomy_database} 2>>log.txt

${cwd}/Utilities/procedure_Orthodb_read_tables.r -l $level -m 1 --odb HOGs_${level}_${taxonomy_database} 2>>log.txt

#taxonomy-constrained tree with RaxML 
if [ $tree == 'raxml' ] 
then
	raxmlHPC-PTHREADS -g $level.phylip.tree -s $level.phylip.data -n $level -m BINCAT -p 33 -T $ncores 2>>log.txt
fi

#order tables by trees
for d in ${ladder[@]}; do
    ${cwd}/Utilities/procedure_Orthodb_order_by_tree.r $level.RData -t $tree -d $d -o Viridiplantae 2>>log.txt &
done
wait


#cotr analysis
for d in ${ladder[@]}; do
    ${cwd}/cotr_transitions.py -m 4 $level.$tree.$d.csv.num 2>>log.txt| ${cwd}/cotr_Fisher.r -pa 1 -p 1e-3 - > $level.$tree.$d.transitions.annotated 2>>log.txt &
done
wait

#results intersections
cat $level.$tree.*.transitions.annotated | perl -anE 'print if ($F[7]>0 and $F[10]<1e-3 or $F[0]=~"Orthogroup")' | cut -f1,2 | sort | uniq -dc| grep -w ${#ladder[@]} | cut -b9- > $level.common.txt 
grep -w -f $level.common.txt $level.$tree.${ladder[0]}.transitions.annotated > intersection.$level.$tree.${ladder[0]}.transitions.annotated

#cluster with mcl
${cwd}/Utilities/procedure_mcl_clusters.pl $level.$tree.${ladder[0]}.csv intersection.$level.$tree.${ladder[0]}.transitions.annotated 2>>log.txt
