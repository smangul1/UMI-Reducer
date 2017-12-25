import sys
import csv
import os
import re
import gzip
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('input', help='<fastq.gz> file')
ap.add_argument('outDir', help='Directory to be created to save <fastq> files. Each <fastq> file is obtained based on the barcode')
optional_arguments = ap.add_argument_group("Optional Inputs")
optional_arguments.add_argument("-min", "--minLength", help="The min length of the tag(read) extcated. Default value is 30 amino acids.", type=int)
optional_arguments.add_argument("-max", "--maxLength", help="The max length of the tag(read) extcated. Default value is 50 amino acids.", type=int)
args = ap.parse_args()





if args.minLength:
    minLength=int(args.minLength)
else:
    minLength=30


if args.maxLength:
    maxLength=int(args.maxLength)
else:
    maxLength=50


if minLength>maxLength:
    print ("ERROR : minLength needs to be smaller than maxLength")
    sys.exit(2)



if not os.path.exists(args.outDir):
    os.makedirs(args.outDir)








def distance(a, b):
    return  sum(map(lambda (x, y): 0 if x == y else 1, zip(a, b)))


#GTGTCAGTCACTTCCAGCGG

def tag_in(read):
    for i in range(14+30,13+50):
        #print i, read[i]
        if distance(read[i:i+20],"GTGTCAGTCACTTCCAGCGG")==1:
            return i
    return -1


#tag_in('CCCATTTGCTATACTAGGACATTTAATGATTGTTGGCAAGAAATGTTCGTGTCAGTCACTTCCAGCGGCAGGAATGCCGAGATCTCTTATGCCGTCTTCTGCTTGAATAAATATCTTCTTTTTTTTCTTTTTTTTTTTTTTTTTTTTTTTC')




#h0 - high quality with ED=0
#h1- high quality with ED=1
#m0 - medium quality with ED=0


#high quality - at least 75% of the bases with Q>30
#medium quality - at least 75% of the bases with Q>20

#NNNNNNNCGATDDD - tag (variable length) - GTGTCAGTCACTTCCAGCGG




fastq_f=args.input
dir_out=args.outDir

statDir=dir_out+"/nReadsStat/"

if not os.path.exists(statDir):
    os.makedirs(statDir)

statFile=statDir+"/nReadsStat.csv"

print (fastq_f)
print (dir_out)

sample_name=os.path.basename(fastq_f).split('.fastq')[0]

print (sample_name)




handle = gzip.open(fastq_f)
seqs={}


barcode=[]
barcode_dict={}

nReads={}


k=0
while True:
        try:
            seq_id = handle.next().strip("\n")
            seq = handle.next().strip("\n")
            t3=handle.next().strip("\n")
            t4=handle.next().strip("\n")
            #print seq
            

            
            if re.search("GTGTCAGTCACTTCCAGCGG", seq):
                p2=re.search("GTGTCAGTCACTTCCAGCGG", seq).start()
                                
                tag=seq[14:p2]
                
                
                #h0 - high quality with ED=0
                k_quality_30=0
                for i in t4[p2:p2+20]:
                    if ord(i)>30:
                        k_quality_30+=1
            
                ratio_30=float(k_quality_30)/20.0
            
                #medium quality - at least 75% of the bases with Q>20
                k_quality_20=0
                for i in t4[p2:p2+20]:
                    if ord(i)>20:
                        k_quality_20+=1
                ratio_20=float(k_quality_20)/20.0
                
                
                
                if len(tag)>=minLength and len(tag)<=maxLength and (ratio_20>=0.75 or ratio_30>=0.75):
                    #17117312_h_0_GTTGCTC_AGCT
                    if seq[7:11] not in barcode:
                        file_barcode=dir_out+sample_name+"_"+seq[7:11]+".fastq"
                        print ("Create",file_barcode)
                        barcode_dict[seq[7:11]] = open(file_barcode, 'w')
                        barcode.append(seq[7:11])
                        nReads[seq[7:11]]=[0,0,0]
                    
                    
                        
                    if ratio_30>=0.75:
                        read_name="@"+str(k)+"_h_0_"+seq[0:7]+"_"+seq[7:11]
                    elif ratio_20>=0.75:
                        read_name="@"+str(k)+"_m_0_"+seq[0:7]+"_"+seq[7:11]
                    k+=1
                    
                
                    barcode_dict[seq[7:11]].write(read_name)
                    barcode_dict[seq[7:11]].write("\n")
                    barcode_dict[seq[7:11]].write(tag)
                    barcode_dict[seq[7:11]].write("\n")
                    barcode_dict[seq[7:11]].write(t3)
                    barcode_dict[seq[7:11]].write("\n")
                    barcode_dict[seq[7:11]].write(t4[14:p2])
                    barcode_dict[seq[7:11]].write("\n")
                    nReads[seq[7:11]][1]+=1
                
            elif tag_in(seq)!=-1:
                p2=tag_in(seq)
                tag=seq[14:p2]
                
                #h0 - high quality with ED=0
                k_quality_30=0
                for i in t4[p2:p2+20]:
                    if ord(i)>30:
                        k_quality_30+=1
                
                ratio_30=float(k_quality_30)/20.0
                
    
                if len(tag)>=minLength and len(tag)<=maxLength and ratio_30>=0.75:
                    #17117312_h_0_GTTGCTC_AGCT
                    if seq[7:11] not in barcode:
                        file_barcode=dir_out+sample_name+"_"+seq[7:11]+".fastq"
                        print ("Create",file_barcode)
                        barcode_dict[seq[7:11]] = open(file_barcode, 'w')
                        barcode.append(seq[7:11])
                        nReads[seq[7:11]]=[0,0,0]
    
                    read_name="@"+str(k)+"_h_1_"+seq[0:7]+"_"+seq[7:11]
                    k+=1
    
    
                    barcode_dict[seq[7:11]].write(read_name)
                    barcode_dict[seq[7:11]].write("\n")
                    barcode_dict[seq[7:11]].write(tag)
                    barcode_dict[seq[7:11]].write("\n")
                    barcode_dict[seq[7:11]].write(t3)
                    barcode_dict[seq[7:11]].write("\n")
                    barcode_dict[seq[7:11]].write(t4[14:p2])
                    barcode_dict[seq[7:11]].write("\n")
                    nReads[seq[7:11]][1]+=1
            else:
                    if seq[7:11] not in barcode:
                        nReads[seq[7:11]]=[0,0,0]
                    nReads[seq[7:11]][2]+=1
                
            
        except StopIteration:
            break
handle.close()

for b in barcode:
    barcode_dict[b].close()

fStat=open(statFile, 'w')
fStat.write("barcode,Total,Good,Bad\n")

for k,v in nReads.items():
    fStat.write(k+","+str(v[1]+v[2])+","+str(v[1])+","+str(v[2]))
    fStat.write("\n")






fStat.close()







