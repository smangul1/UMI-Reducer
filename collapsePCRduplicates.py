import pysam
import sys
import collections
import numpy as np
import random
import csv
import argparse
import matplotlib
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt


#updated 09/04/2015


def editDistance(s1,s2):
    k=0
    for i in range(0,len(s1)):
        if s1[i]!=s2[i]:
            k+=1
    return k


ap = argparse.ArgumentParser()
ap.add_argument('inbam', help='Mapped reads in bam format')
ap.add_argument('outbam', help='Output file to save reads after collapsing PCR duplicates')
#ap.add_argument('--testN', type=int,
#                help='Run a test using only the first N features, and then '
#                'print out some example feature IDs and their attributes')
ap.add_argument('--m', action='store_true',help='Save multi-mapped reads')

#cmd https://gist.github.com/daler/ec481811a44b3aa469f3

args = ap.parse_args()





bam=args.inbam
out=args.outbam



chr_list=[]

for i in range(1,20):
    chr_list.append(str(i))

chr_list.append('X')
chr_list.append('Y')
chr_list.append('MT')

print chr

position=[]
position_all_uniq=[]





samfile = pysam.AlignmentFile(bam, "rb" )


dict= {}


mappedReads=[]
numberReadsUnique=0

numberReadsUniquePlusMultiMapped=0

numberReadsUnique_covGreated1=0
numberReadsUnique_filtered=0
readLength=[]
readLength_filtered=[]
before=0
after=0







bam_header = pysam.Samfile(bam, 'rb').header



outfile = pysam.AlignmentFile(out, "wb", header=bam_header)








for chr in chr_list:
    dict.clear()
    position[:]=[]
    print "----------chr",chr
    for read in samfile.fetch(chr):
        mappedReads.append(read.query_name)
        
        
        if args.m:
            if read.mapq==50:
                numberReadsUnique+=1
            numberReadsUniquePlusMultiMapped+=1
            position.append(read.reference_start)
            readLength.append(len(read.query_sequence))
        else:
            if read.mapq==50:
                numberReadsUnique+=1
                position.append(read.reference_start)
                readLength.append(len(read.query_sequence))


    print "numberReadsUnique",numberReadsUnique
    print "numberReadsUniquePlusMultiMapped",numberReadsUniquePlusMultiMapped



    counter_chr=collections.Counter(position)
    position_all_uniq+=position
    print "Number of position with #reads staring >=1", len(position)

    
    for key,val in counter_chr.items():
        #print key,val

        if val==1:
            for read in samfile.fetch(chr,key,key+1):
                if read.reference_start==key:
                    outfile.write(read)
                    readLength_filtered.append(len(read.query_sequence))
                    numberReadsUnique_filtered+=1
     
    



        elif val>1:
            Read=[]
            Read[:]=[]
            for read in samfile.fetch(chr,key,key+1):
                if read.reference_start==key:
                    Read.append(read)
            
        #print key,val
           
            
            
            exclude_list=[]
            exclude_list[:]=[]
            
            for i in range(0,val):
                for j in range(0,val):
                    if i<j:
                        #print i,j
                        #print Read[dict[key]][i].query_name, Read[dict[key]][i].query_sequence
                        #print Read[dict[key]][j].query_name, Read[dict[key]][j].query_sequence

                        #print Read[dict[key]][0].query_name
                        #print Read[dict[key]][0].query_name.split("_")[3]+"_"+Read[dict[key]][0].query_sequence
                        if Read[i].query_name.split("_")[3]+"_"+Read[i].query_sequence==Read[j].query_name.split("_")[3]+"_"+Read[j].query_sequence:
                            #print i,j
                            exclude_list.append(j)
                            #print "EQUAL"
                            #print Read[dict[key]][i].query_name, Read[dict[key]][i].query_sequence
                            #print Read[dict[key]][j].query_name, Read[dict[key]][j].query_sequence
                            #print "to exlude",j
        
        
            for i in range(0,val):
                numberReadsUnique_covGreated1+=1
                if i not in exclude_list:
                    #print "to include",i
                    readLength_filtered.append(len(Read[i].query_sequence))
                    numberReadsUnique_filtered+=1
                    outfile.write(Read[i])
                    

outfile.close()



#-----------------------
#statistics


header=[]

header.append('sample')
header.append('Number of mapped reads')
header.append('Number of reads after removing reads mapped to multiple locations in the genome')
header.append('Number of reads after collapsing PCR dublicated')

nr=[]

nr.append(out.split('.')[0])
nr.append(len(set(mappedReads)))
nr.append(numberReadsUnique)

nr.append(numberReadsUniquePlusMultiMapped)

nr.append(numberReadsUnique_filtered)




stat_f=out.split('.')[0]+'.number_of_reads_stat'

with open(stat_f, 'w') as fp:
    a = csv.writer(fp, delimiter=',')
    a.writerow(header)
    a.writerow(nr)




counter=collections.Counter(position_all_uniq)
position_all_uniq=set(position_all_uniq)
print "Number of position with #reads staring >=1", len(position_all_uniq)





#-----------------------
#save as a histogram - counter_length_filtered

x1=[]
xbins1=[]

counter_length_filtered=collections.Counter(readLength_filtered)
for key,val in counter_length_filtered.items():
    #print key,val
    xbins1.append(val)
    x1.append(key)



plot1=out.split('.')[0]+'.readLengthPCRDuplicatesCollapsed.png'
print "save to",plot1



plt.title('Length of reads after collapsing PCR duplicates')

plt.bar(x1,xbins1)
plt.savefig(plot1)



#save as a histogram - counter_length_uniq
x2=[]
xbins2=[]

counter_length=collections.Counter(readLength)
for key,val in counter_length.items():
    #print key,val
    xbins2.append(val)
    x2.append(key)

plot2=out.split('.')[0]+'.readLengthMultiMappedReadsRemoved.png'
print "save to",plot2



plt.title('Length of after removing reads mapped to multiple locations in the genome')
plt.bar(x2,xbins2)
plt.savefig(plot2)





print "DONE!"


