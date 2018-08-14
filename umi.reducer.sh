#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1

parser.add_argument('bam',help='bam with mapped and unmapped reads, BAM file needs to be indexed')
#parser.add_argument('-chr', '--chr', default='', type=str,help='Use this option to specify the format of chromosome name in the bam file. By default it is : 1,2,3,...,X. In case the format is : chr1, chr2,..,chrX, use : -chr chr)')
EOF


DIR_CODE=`dirname $(readlink -f "$0")`

if [ -d "$DIR_CODE/tools/MiniConda/bin" ]; then
    echo "Add MiniConda to PATH if it's available"
    export PATH="$DIR_CODE/tools/MiniConda/bin:$PATH"
fi


#filename=$(basename "$BAM")
extension="${BAM##*.}"
filename="${BAM%.*}"
prefix=${filename##*/}



dirSource="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Start analysis"


#collapsePCRduplicates.py
python ${dirSource}/collapsePCRduplicates.py  ${BAM} ${filename}_PCRduplicates_removed.bam
samtools sort ${filename}_PCRduplicates_removed.bam >${filename}_PCRduplicates_removed.sort.bam

rm ${filename}_PCRduplicates_removed.bam
samtools index ${filename}_PCRduplicates_removed.sort.bam

echo "Run genomic.categories.py ..."


python ${dirSource}/mRNA.abundance.py --perCategory --mouse --multi ${filename}_PCRduplicates_removed.sort.bam ${filename}









#echo ${SPECIES[@]}
#echo required infile: $BAM


echo "DONE!"
