

if [ $# -lt 1 ]
then


echo "[1] - Directory with fastq files"
#echo "[3] - Result directory. Bam file with collapsed PCR dublicates are saved here"
exit 1
fi

dirSource="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $dirSource


dir=$1

cd $dir
ls *fastq | awk -F ".fastq" '{print $1}' >samples.txt


while read line
do
mkdir $line
cd $line

echo ". /u/local/Modules/default/init/modules.sh">run_${line}.sh
echo "module load python/2.7.3">>run_${line}.sh

echo "${dirSource}/run.hisat2.tuned.SE.mouse.sh ${dir}/${line}.fastq ${dir}/${line}/${line}/">>run_${line}.sh
echo "${dirSource}/umi-reducer.sh ${dir}/${line}..sort.bam">>run_${line}.sh


cd ..

done<samples.txt

exit 1

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






