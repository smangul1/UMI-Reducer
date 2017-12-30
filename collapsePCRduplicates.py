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
ap.add_argument("--MAPQ",help="mapping quality to extract uniquely-mapped reads. Only reads with this MAPQ will be selected. By default all reads will be considered", default=60,type=int)
ap.add_argument("--u",help="a binary flag used to indicate that only uniquely mapped reads will be considered. By default uniquely mapped reads are defined as reads with MAPQ=60",action="store_true")
ap.add_argument("--chr",help="Number of autosomes. By default 20", default=20,type=int)
ap.add_argument("--end",help="a binary flag used to indicate that end position is considered. This is beta version use with caution",action="store_true")
ap.add_argument('--e', action='store_true',help='UMI-tools example format')



args = ap.parse_args()








bam=args.inbam
out=args.outbam



chr_list=[]

if args.e:
    for i in range(1,20):
        chr_list.append('chr'+str(i))

    chr_list.append('chrX')
    chr_list.append('chrY')
    chr_list.append('chrM')

else:
    for i in range(1,args.chr):
        chr_list.append(str(i))

    chr_list.append('X')
    chr_list.append('Y')
    chr_list.append('MT')

print (chr)

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




readSet=set()


bam_header = pysam.Samfile(bam, 'rb').header



outfile = pysam.AlignmentFile(out, "wb", header=bam_header)



print ("Open ",bam, "via pysam")




for chr in chr_list:
    dict.clear()
    position[:]=[]
    if args.e:
        print "----------",chr
    else:
        print ("----------chr",chr)
    for read in samfile.fetch(chr):
        mappedReads.append(read.query_name)

        if args.u:
            if read.mapq==args.MAPQ:
                numberReadsUnique += 1
                position.append(read.reference_start)
                readLength.append(len(read.query_sequence))
        else:
            numberReadsUniquePlusMultiMapped+=1
            position.append(read.reference_start)
            readLength.append(len(str(read.query_sequence)))



    print ("numberReadsUnique",numberReadsUnique)
    print ("numberReadsUniquePlusMultiMapped",numberReadsUniquePlusMultiMapped)



    counter_chr=collections.Counter(position)
    position_all_uniq+=position
    print ("Number of position with #reads staring >=1", len(position))

    count=0

    print  ("Processing", len(counter_chr.items()), "items")

    for key,val in counter_chr.items():
        if count%10000==1:
            print (count)
        count+=1
        
        if val==1:

            for read in samfile.fetch(chr,key,key+1):
                if read.reference_start==key:

                    outfile.write(read)
                    readLength_filtered.append(len(str(read.query_sequence)))
                    numberReadsUnique_filtered+=1
                    readSet.add(read.query_name)
     
    




        elif val>1:
            setReads=set()
            setReads.clear()
    
            if val>1000:
                print (val, chr,key)
            Read=[]
            Read[:]=[]
            for read in samfile.fetch(chr,key,key+1):
                if read.reference_start==key:
                    Read.append(read)
                    if args.end:
                        #UMI-tools format
                        if args.e:
                            setReads.add(read.query_name.split("_")[1] + "_" + str(read.query_sequence))
                        #default format
                        else:
                            setReads.add(read.query_name.split("_")[3] + "_" + read.query_sequence)
                    else:
                        #UMI-tools format
                        if args.e:
                            setReads.add(read.query_name.split("_")[1]+"_"+str(read.query_sequence))
                        #default format
                        else:
                            setReads.add(read.query_name.split("_")[3]+"_"+read.query_sequence+"_"+str(read.reference_end))


                    #print (read.query_sequence,read.query_name,read.query_name.split("_")[3],read.reference_start,read.reference_end)

            

            notsetReads=set()
            notsetReads.clear()
            numberReadsUnique_covGreated1+=len(setReads)
            for i in range(0,val):
                if args.end:
                    #UMI-tools format
                    if args.e:
                        extended_read_name = Read[i].query_name.split("_")[1] + "_" + str(Read[i].query_sequence)+"_"+str(Read[i].reference_end)
                    #default format
                    else:
                        extended_read_name = Read[i].query_name.split("_")[3] + "_" + Read[i].query_sequence+"_"+str(Read[i].reference_end)
                else:
                    #UMI-tools format
                    if args.e:
                        extended_read_name=Read[i].query_name.split("_")[1]+"_"+Read[i].query_sequence
                    #default format
                    else:
                        extended_read_name=Read[i].query_name.split("_")[3]+"_"+Read[i].query_sequence

                if extended_read_name in setReads and extended_read_name not in notsetReads:
                        outfile.write(Read[i])
                        numberReadsUnique_filtered+=1
                        readLength_filtered.append(len(str(Read[i].query_sequence)))


                        if args.end:
                            #UMI-tools format
                            if args.e:
                                notsetReads.add(Read[i].query_name.split("_")[1] + "_" + Read[i].query_sequence+"_"+str(Read[i].reference_end))
                            #default format
                            else:
                                notsetReads.add(Read[i].query_name.split("_")[3] + "_" + Read[i].query_sequence+"_"+str(Read[i].reference_end))
                        else:
                            #UMI-tools format
                            if args.e:
                                notsetReads.add(Read[i].query_name.split("_")[1]+"_"+Read[i].query_sequence)
                            #default format
                            else:
                                notsetReads.add(Read[i].query_name.split("_")[3]+"_"+Read[i].query_sequence)

                        readSet.add(Read[i].query_name)



                    

outfile.close()



#-----------------------
#statistics


header=[]

header.append('sample')
header.append('Number of mapped reads')
header.append('Number of reads mapped to unique location (UNIQUE reads)')
header.append('Number of reads after collapsing PCR dublicated (each read is present once) ')


nr=[]

nr.append(out.split('.')[0])
nr.append(len(set(mappedReads)))
nr.append(numberReadsUnique)
nr.append(len(readSet))




stat_f=out.split('.')[0]+'.number_of_reads_stat'
print ("Save to ", stat_f)

with open(stat_f, 'w') as fp:
    a = csv.writer(fp, delimiter=',')
    a.writerow(header)
    a.writerow(nr)




counter=collections.Counter(position_all_uniq)
position_all_uniq=set(position_all_uniq)
print ("Number of position with #reads staring >=1", len(position_all_uniq))




#-----------------------
#save as a histogram - counter_length_filtered

x1=[]
xbins1=[]

counter_length_filtered=collections.Counter(readLength_filtered)
for key,val in counter_length_filtered.items():
    xbins1.append(val)
    x1.append(key)



plot1=out.split('.')[0]+'.readLengthPCRDuplicatesCollapsed.png'
print ("save to",plot1)



plt.title('Length of reads after collapsing PCR duplicates')

plt.bar(x1,xbins1)
plt.savefig(plot1)



#save as a histogram - counter_length_uniq
x2=[]
xbins2=[]

counter_length=collections.Counter(readLength)
for key,val in counter_length.items():
    xbins2.append(val)
    x2.append(key)

plot2=out.split('.')[0]+'.readLengthBeforePCRduplicates.png'
print ("save to",plot2)



plt.title('Length before collapsing PCR duplicates')
plt.bar(x2,xbins2)
plt.savefig(plot2)


samfile.close()


print ("DONE!")


