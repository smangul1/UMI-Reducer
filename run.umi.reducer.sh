

if [ $# -lt 1 ]
then


echo "[1] - Directory with fastq files"
exit 1
fi

dirSource="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo $dirSource



dir=$1
cd $dir
ls *fastq | awk -F ".fastq" '{print $1}' >samples.txt





while read line
do
#mkdir $line
#cd $line

echo ". /u/local/Modules/default/init/modules.sh">run_${line}.sh
echo "module load python/2.7.3">>run_${line}.sh

echo "${dirSource}/run.hisat2.tuned.SE.sh ${dir}/${line}.fastq ${dir}/${line}/">>run_${line}.sh
echo "${dirSource}/umi.reducer.sh ${dir}/${line}/${line}.sort.bam">>run_${line}.sh




done<samples.txt



cd $dir
ls *sh | awk '{i+=1;print "if [ $1 == "i" ];then ./"$1" ;fi"}' > myFunc.sh
cp ${dirSource}/myFuncFastWrapper.sh  ./

chmod 755 *sh
chmod 755 */*sh

n=$(cat myFunc.sh | wc -l)
cd $dir
qsub -cwd -V -N umi.reducer -l h_data=16G,time=05:00:00 -t 1-${n}:1 myFuncFastWrapper.sh






