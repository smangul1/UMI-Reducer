

if [ $# -lt 3 ]
then

echo "[1] - File with sample names"
echo "[2] - Directory with bam files"
echo "[3] - Result directory. Bam file with collapsed PCR dublicates are saved here"
exit 1
fi

dirSource="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $dirSource



dir=$3
mkdir $dir
cd $dir

while read line
do
mkdir $line
cd $line

echo ". /u/local/Modules/default/init/modules.sh">run_${line}.sh
echo "module load python/2.7.3">>run_${line}.sh
echo "module load samtools">>run_${line}.sh
echo "samtools index ${2}/${line}.bam">>run_${line}.sh

echo "python ${dirSource}/collapsePCRduplicates.py --m ${2}/${line}.bam ${dir}/${line}/${line}_PCRDuplicatesCollapsed.bam">>run_${line}.sh
echo ". /u/local/Modules/default/init/modules.sh">>run_${line}.sh
echo "module load samtools">>run_${line}.sh
echo "samtools sort ${dir}/${line}/${line}_PCRDuplicatesCollapsed.bam -o ${dir}/${line}/${line}_PCRDuplicatesCollapsed_sort.bam">>run_${line}.sh
echo "samtools index ${dir}/${line}/${line}_PCRDuplicatesCollapsed_sort.bam">>run_${line}.sh
echo "rm -rf ${dir}/${line}/${line}_PCRDuplicatesCollapsed.bam">>run_${line}.sh 



cd ..
done<$1


cd $dir
ls */*sh | awk '{i+=1;print "if [ $1 == "i" ];then ./"$1" ;fi"}' > myFunc.sh
cp ${dirSource}/myFuncFastWrapper.sh  ./

chmod 755 *sh
chmod 755 */*sh

n=$(cat myFunc.sh | wc -l)
echo "qsub -cwd -V -N sarahSeq -l h_data=16G,express,time=24:00:00 -t 1-${n}:1 myFuncFastWrapper.sh
"

cd $dir
qsub -cwd -V -N sarahSeq -l h_data=16G,express,time=24:00:00 -t 1-${n}:1 myFuncFastWrapper.sh






