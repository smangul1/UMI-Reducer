#!/bin/bash

if [ $# -lt 1 ]
then
#echo "[1] - file with samples located in the "
echo "[1] dir with bams"
echo "[2] outDir"
exit 1
fi

gprofileOut=$2
mkdir $gprofileOut

cwd=$(pwd)

cd $1
ls *bam | awk -F ".bam" '{print $1}' >${cwd}/samples4gprofile.txt

cd $cwd

while read line ;do echo ". /u/local/Modules/default/init/modules.sh" >run_gprofile_${line}.sh; echo "module load python/2.7.3" >>run_gprofile_${line}.sh; echo "python /u/home/s/serghei/collab/code/rop/gprofile.py ${1}/${line}.bam ${gprofileOut}/${line}.csv --perCategory --mouse">>run_gprofile_${line}.sh;done<samples4gprofile.txt

ls run_gprofile_*sh | awk '{i+=1;print "qsub -cwd -V -N gprofile"i" -l h_data=16G,time=24:00:00 "$1}' > all.sh

chmod u+x run_gprofile_*sh
chmod u+x all.sh
nohup ./all.sh &

echo "Jobs have been submitted!!!"




