#!/bin/bash

source $(dirname $0)/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1

parser.add_argument('bam',help='bam with mapped and unmapped reads, BAM file needs to be indexed')
#parser.add_argument('-chr', '--chr', default='', type=str,help='Use this option to specify the format of chromosome name in the bam file. By default it is : 1,2,3,...,X. In case the format is : chr1, chr2,..,chrX, use : -chr chr)')
EOF

#filename=$(basename "$BAM")
extension="${BAM##*.}"
filename="${BAM%.*}"
echo $filename


python collapsePCRduplicates.py --m ${BAM} ${filename}_PCRduplicates_removed.bam
tools/samtools-1.3/samtools sort ${filename}_PCRduplicates_removed.bam >${filename}_PCRduplicates_removed.sort.bam

rm ${filename}_PCRduplicates_removed.bam
tools/samtools-1.3/samtools index ${filename}_PCRduplicates_removed.sort.bam

 




#echo ${SPECIES[@]}
#echo required infile: $BAM
