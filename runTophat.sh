#!/bin/bash




if [ $# -lt 6 ]
then
echo "********************************************************************"
echo "RNA-Seq data analyses on hoffman2 : tophat "
echo "For questions or suggestions contact:"
echo "Serghei Mangul smangul@ucla.edu"
echo "********************************************************************"
echo ""

echo ""
echo "1 <fastq_dir>  - dir with  fastq files"
echo "2 <wdir> - working directory, directory where results from rnaseq anlyses will be saved, for each samples separately directory will be created"
echo "3 <bowtie2_index - location of bowtie2 index, e.g. for human /u/home/eeskin/serghei/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome for Drosophila_melanogaster ~/scratch/Drosophila_melanogaster/UCSC/dm3/Sequence/Bowtie2Index/genome"
echo "4 tophat options if more then 2 arguments in the \" :::::: e.g. \"-G ~/scratch/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf\""
echo "5 - file with the fastq names"
echo "6 - dir whre all the scripts of sarahSeq are stored"

echo "--------------------------------------"
exit 1
fi


##Sample_12R0008_1.fastq



dir=$2
bowtie2_index=$3
mkdir $2
cd $dir

echo "tophat options"

tophat_options=$4
echo $tophat_options

samples=$5




########1M reads
p=4000000
#p=400

echo "Advanced options for tophat ..."
echo $tophat_options





while read sdir
do

echo $sdir

ln -s ${1}/${sdir}.fastq ./


f1=${dir}/${sdir}.fastq

rm -rf $sdir
mkdir $sdir
cd $sdir
echo "cd $sdir"
pwd


echo $f1

n=$(cat $f1 | wc -l)
echo $n > nr.txt
m=$(echo $n| awk  '{ rounded = sprintf("%.0f", $1/'${p}'); print rounded }') #million
echo $m >m.txt
m1=$(echo $m | awk '{print $1-1}')
echo $m1

if [ $m == 0 ]
then
i=0
pwd

echo "#!/bin/bash">${i}.sh
echo "cd ${dir}/${sdir}">>${i}.sh
#echo "cat ${f1} | awk '{if((k>=(0+$i)*$p) && (k<(1+$i)*$p)) print;k++ }' > ${dir}/${sdir}/${sdir}_1_${i}.fastq" >>${i}.sh
#-p 8
echo ". /u/local/Modules/default/init/modules.sh" >>${i}.sh
echo "module load tophat/2.0.9" >>${i}.sh
echo "module load bowtie2/2.1.0" >>${i}.sh
echo "module load samtools"  >>${i}.sh
echo "tophat -v" >>${i}.sh
echo "bowtie -v" >>${i}.sh
echo "tophat -o tophat_out_${i} $tophat_options ${bowtie2_index} ${dir}/${sdir}.fastq" >>${i}.sh
echo "mv ${dir}/${sdir}/tophat_out_${i}/accepted_hits.bam ${dir}/${sdir}/${i}_mapped.bam" >>${i}.sh
echo "mv ${dir}/${sdir}/tophat_out_${i}/unmapped.bam ${dir}/${sdir}/${i}_unmapped.bam" >>${i}.sh
ls
echo "0" >m.txt
cd ..
continue

fi


for((i=0;i<$((m-1));i++))
do

echo "Split reads"
echo $i

echo "#!/bin/bash">${i}.sh
echo "cd ${dir}/${sdir}">>${i}.sh
echo "cat ${f1} | awk '{if((k>=(0+$i)*$p) && (k<(1+$i)*$p)) print;k++ }' > ${dir}/${sdir}/${sdir}_1_${i}.fastq" >>${i}.sh
#-p 8
echo ". /u/local/Modules/default/init/modules.sh" >>${i}.sh
echo "module load tophat/2.0.9" >>${i}.sh
echo "module load bowtie2/2.1.0" >>${i}.sh
echo "module load samtools"  >>${i}.sh
echo "tophat -v" >>${i}.sh
echo "bowtie -v" >>${i}.sh
echo "tophat -o tophat_out_${i} $tophat_options ${bowtie2_index} ${dir}/${sdir}/${sdir}_1_${i}.fastq" >>${i}.sh
echo "mv ${dir}/${sdir}/tophat_out_${i}/accepted_hits.bam ${dir}/${sdir}/${i}_mapped.bam" >>${i}.sh
echo "mv ${dir}/${sdir}/tophat_out_${i}/unmapped.bam ${dir}/${sdir}/${i}_unmapped.bam" >>${i}.sh
#t    echo "rm ${sdir}_1_${i}.fastq" >>${i}.sh
#t    echo "rm ${sdir}_2_${i}.fastq" >>${i}.sh
#t    echo "rm -rf tophat_out_${i}/logs" >>${i}.sh
#t    echo "rm -rf tophat_out_${i}/tmp" >>${i}.sh




done

echo "#!/bin/bash">${m1}.sh

echo "cd ${dir}/${sdir}">>${m1}.sh
echo "cat ${f1} | awk '{if(k>=($m-1)*$p) print;k++ }' > ${dir}/${sdir}/${sdir}_1_${m1}.fastq" >>${m1}.sh
#-p 8
echo ". /u/local/Modules/default/init/modules.sh" >>${m1}.sh
echo "module load tophat/2.0.9" >>${m1}.sh
echo "module load bowtie2/2.1.0" >>${m1}.sh
echo "module load samtools" >>${m1}.sh
echo "tophat -v" >>${m1}.sh
echo "bowtie -v" >>${m1}.sh
echo "tophat  -o tophat_out_${m1}  $tophat_options  ${bowtie2_index} ${dir}/${sdir}/${sdir}_1_${m1}.fastq" >>${m1}.sh
echo "mv ${dir}/${sdir}/tophat_out_${m1}/accepted_hits.bam ${dir}/${sdir}/${m1}_mapped.bam" >>${m1}.sh
echo "mv ${dir}/${sdir}/tophat_out_${m1}/unmapped.bam ${dir}/${sdir}/${m1}_unmapped.bam" >>${m1}.sh

#t echo "rm ${sdir}_1_${m1}.fastq">>${m1}.sh
#t echo "rm ${sdir}_2_${m1}.fastq" >>${m1}.sh
#t echo "rm -rf tophat_out_${m1}/logs" >>${m1}.sh
#t echo "rm -rf tophat_out_${m1}/tmp" >>${m1}.sh






cd ..
pwd


done<$samples

cp $samples samples.txt

echo "================"
pwd

#chmod 755 myFunc.sh
#chmod 755 *sh
#chmod 755 */*sh

ls */*sh | awk '{i+=1;print "if [ $1 == "i" ];then ./"$1" ;fi"}' > myFunc.sh
cp ${6}//myFuncFastWrapper.sh ./

chmod 755 myFunc.sh
chmod 755 *sh
chmod 755 */*sh

n=$(cat myFunc.sh | wc -l)
echo "qsub -cwd -V -N tophatA -l h_data=8G,express,time=10:00:00 -t 1-${n}:1 myFuncFastWrapper.sh"

qsub -cwd -V -N tophatA -l h_data=8G,express,time=10:00:00 -t 1-${n}:1 myFuncFastWrapper.sh

echo "The jobs were submitted!!!!"

#echo "Please ak 8G from qsub
echo "done!




