#!/bin/bash

AUTHOR="Serghei Mangul"



################################################################
##########          The main template script          ##########
################################################################

toolName="hisat2.tuned"
toolPath="/u/home/s/svandrie/code/import/hisat2-2.1.0/hisat2"
index=/u/home/s/svandrie/project-kmartin/Mus_musculus/NCBI/build37.2/Sequence/hisat2/Mus.musculus.build37.2.hisat2.index


#index="/u/home/h/harryyan/project-eeskin/utilities/hisat2-2.1.0/ref_genome/grch38/genome"






if [ $# -lt 2 ]
    then
    echo "********************************************************************"
    echo "This script was written by Serghei Mangul"
    echo "********************************************************************"
    echo ""
    echo "1 <input>   - R.fastq"
    #echo "2 <input2>   - R2.fastq"
    echo "2 <outdir>  - dir to save the output"
    echo "--------------------------------------"
    exit 1
    fi


DIR_CODE=`dirname $(readlink -f "$0")`

if [ -d "$DIR_CODE/tools/MiniConda/bin" ]; then
    echo "Add MiniConda to PATH if it's available"
    export PATH="$DIR_CODE/tools/MiniConda/bin:$PATH"
fi


# mandatory part
input1=$1
#input2=$2
outdir=$2
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


#Convert to absolute paths.
input1=`readlink -m "$input1"`
outdir=`readlink -m "$outdir"`

# STEP 0 - create output directory if it does not exist

mkdir $outdir
pwd=$PWD
cd $outdir
outdir=$PWD
cd $pwd
logfile=$outdir/report_$(basename ${input1%.*})_${toolName}.log


echo $logfile


# -----------------------------------------------------

echo "START" >> $logfile
# STEP 1 - prepare input if necessary (ATTENTION: TOOL SPECIFIC PART!)
# -----------------------------------


# STEP 2 - run the tool (ATTENTION: TOOL SPECIFIC PART!)

now="$(date)"
printf "%s --- RUNNING %s\n" "$now" $toolName >> $logfile

# run the command
res1=$(date +%s.%N)


. /u/local/Modules/default/init/modules.sh
#module load samtools


hisat2 -x $index -U $input1 --end-to-end -N 1 -L 20 -i S,1,0.5 -D 25 -R 5 --pen-noncansplice 12 --mp 1,0 --sp 3,0 --time --reorder | samtools view -F 4 -bS - >$outdir/$(basename ${input1%.*}).bam 2>>$logfile


samtools sort $outdir/$(basename ${input1%.*}).bam >$outdir/$(basename ${input1%.*}).sort.bam
rm -fr $outdir/${toolName}_$(basename ${input1%.*}).bam
samtools index $outdir/$(basename ${input1%.*}).sort.bam

#tools/samtools-1.3/samtools view -f 0x4 -bh $outdir/${toolName}_$(basename ${input1%.*}).bam | samtools bam2fq - >$outdir/${toolName}_$(basename ${input1%.*})_unmapped.fastq 2>>$logfile

#hisat2 --threads 16 --end-to-end -N <NUM_MISMATCH> -L <SEED_LENGTH> -i S,1,<SEED_INTERVAL> -D <SEED_EXTENSION> -R <RE_SEED> --pen-noncansplice <PENALITY_NONCANONICAL> --mp <MAX_MISMATCH_PENALITY>,<MIN_MISMATCH_PENALITY> --sp <MAX_SOFTCLIPPING_PENALITY>,<MIN_SOFTCLIPPING_PENALITY>--time --reorder --known-splicesite-infile <output index path>/<genome name>.splicesites.txt --novel-splicesite-outfile splicesites.novel.txt --novel-splicesite-infile splicesites.novel.txt -f -x <index name> -1 <read file 1> -2 <read file 2> -S <output sam file>
# default 1 20 0.5 25 5 12 1 0 3 0



res2=$(date +%s.%N)
dt=$(echo "$res2 - $res1" | bc)
dd=$(echo "$dt/86400" | bc)
dt2=$(echo "$dt-86400*$dd" | bc)
dh=$(echo "$dt2/3600" | bc)
dt3=$(echo "$dt2-3600*$dh" | bc)
dm=$(echo "$dt3/60" | bc)
ds=$(echo "$dt3-60*$dm" | bc)
now="$(date)"
printf "%s --- TOTAL RUNTIME: %d:%02d:%02d:%02.4f\n" "$now" $dd $dh $dm $ds >> $logfile

now="$(date)"
printf "%s --- FINISHED RUNNING %s %s\n" "$now" $toolName >> $logfile

# ---------------------




# STEP 3 - transform output if necessary (ATTENTION: TOOL SPECIFIC PART!)



now="$(date)"
printf "%s --- TRANSFORMING OUTPUT\n" "$now" >> $logfile


#cat $outdir/one_output_file.fastq | gzip > $outdir/${toolName}_$(basename ${input%.*})_${kmer}.corrected.fastq.gz

now="$(date)"
printf "%s --- TRANSFORMING OUTPUT DONE\n" "$now" >> $logfile

# remove intermediate files
#rm $outdir/one_output_file.fastq


# --------------------------------------



printf "DONE" >> $logfile
