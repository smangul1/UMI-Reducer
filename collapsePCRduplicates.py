import pysam
import sys
import collections
import numpy as np
import random
import csv
import argparse
#import matplotlib
#matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
#import matplotlib.pyplot as plt
import os

#updated 09/04/2015





def editDistance(s1,s2):
    k=0
    for i in range(0,len(s1)):
        if s1[i]!=s2[i]:
            k+=1
    return k

def considerPos():
	if args.end:
		pos = str(read.reference_start)
	else:
		pos = str(read.cigarstring)
	return pos

def considerPos2():
	if args.end:
		pos = str(Read[i].reference_start)
	else:
		pos = str(Read[i].cigarstring)
	return pos

ap = argparse.ArgumentParser()
ap.add_argument('inbam', help='Mapped reads in bam format')
ap.add_argument('outbam', help='Output file to save reads after collapsing PCR duplicates')
ap.add_argument("--MAPQ",help="mapping quality to extract uniquely-mapped reads. Only reads with this MAPQ will be selected. By default all reads will be considered", default=60,type=int)
ap.add_argument("--u",help="a binary flag used to indicate that only uniquely mapped reads will be considered. By default uniquely mapped reads are defined as reads with MAPQ=60",action="store_true")
ap.add_argument("--chr",help="Number of autosomes. By default 20", default=20,type=int)
ap.add_argument('--e', action='store_true',help='UMI-tools example format')
ap.add_argument('--end',action='store_true',help='End position of read will not be considered')


args = ap.parse_args()








bam=args.inbam
out=args.outbam
path=os.path.dirname(out)
base=os.path.basename(out)
prefix=os.path.splitext(base)[0]

chr_list=[]





position=[]
position_all_uniq=[]






samfile = pysam.Samfile(bam, "rb" )

#extract chr names frm BAm file
chr_list=[]
for i in samfile.header['SQ']:
    chr_list.append(i['SN'])


print ("List of chromosomes extracted from BAM")
print (chr_list)

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



outfile = pysam.Samfile(out, "wb", header=bam_header)



print ("Open ",bam, "using pysam")



numberReadsUniqueGlobal=0 # across all chromosomes

for chr in chr_list:

    numberReadsUnique=0
    dict.clear()
    position[:]=[]




    for read in samfile.fetch(chr):
        mappedReads.append(read.qname)

        if read.get_tag("NH") == 1:
            numberReadsUnique += 1
            numberReadsUniqueGlobal+=1

        if args.u:
            if read.get_tag("NH")==1:
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
    print ("Number of position with #reads staring >=1", len(set(position)))

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
    


            Read=[]
            Read[:]=[]
            for read in samfile.fetch(chr,key,key+1):
                if read.reference_start==key:
                    Read.append(read)

                    #UMI-tools format
                    if args.e:
                        setReads.add(read.query_name.split("_")[1]+"_"+considerPos())
                    #default format
                    else:
                        setReads.add(read.query_name.split("_")[3]+"_"+considerPos())




            #for example this is setReads: set(['ACA_31M_168', 'AAA_31M_168'])

            notsetReads=set()
            notsetReads.clear()
            numberReadsUnique_covGreated1+=len(setReads)
            for i in range(0,val):

                #UMI-tools format
                if args.e:
                    extended_read_name=Read[i].query_name.split("_")[1]+"_"+considerPos2()
                    #default format
                else:
                    extended_read_name=Read[i].query_name.split("_")[3]+"_"+considerPos2()



                if extended_read_name in setReads and extended_read_name not in notsetReads:
                        outfile.write(Read[i])
                        numberReadsUnique_filtered+=1
                        readLength_filtered.append(len(str(Read[i].query_sequence)))



                        #UMI-tools format
                        if args.e:
                            notsetReads.add(Read[i].query_name.split("_")[1]+"_"+considerPos2())
                        #default format
                        else:
                            notsetReads.add(Read[i].query_name.split("_")[3]+"_"+considerPos2())



                        readSet.add(Read[i].query_name)





outfile.close()



#-----------------------
#statistics


print ('Number of mapped reads',len(set(mappedReads)))
print ('Number of reads mapped to unique location (UNIQUE reads)',numberReadsUniqueGlobal)
print ('Number of reads after collapsing PCR dublicated (each read is present once)',len(readSet))


header=[]

header.append('sample')
header.append('Number of mapped reads')
header.append('Number of reads mapped to unique location (UNIQUE reads)')
header.append('Number of reads after collapsing PCR dublicated (each read is present once) ')


nr=[]

nr.append(prefix.replace('_PCRduplicates_removed',''))
nr.append(len(set(mappedReads)))
nr.append(numberReadsUniqueGlobal)
nr.append(len(readSet))




stat_f=path+"/"+prefix+'.number_of_reads_stat'
print ("Save to ", stat_f)

#with open(stat_f, 'w') as fp:
   # a = csv.writer(fp, delimiter=',')
   # a.writerow(header)
   # a.writerow(nr)




samfile.close()


print "Create a new BAM file", args.outbam


counter=collections.Counter(position_all_uniq)
position_all_uniq=set(position_all_uniq)
print ("Number of position with #reads staring >=1", len(position_all_uniq))



plot1=path+"/"+prefix+'.readLength.csv'
print ("Save to",plot1)



f=open(plot1, "w")
f.write("flag,l"+"\n")
    
for l in readLength_filtered:
    f.write("PCRDuplicatesCollapsed,"+str(l) +"\n")

f.close()


f=open(plot1, "a")
for l in readLength:
    f.write("BeforePCRduplicates,"+str(l) +"\n")







print ("DONE!")


