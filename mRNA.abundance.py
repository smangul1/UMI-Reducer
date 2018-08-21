import pysam
import sys
import csv
import os
import argparse
from numpy.random import choice
import random
from quicksect import IntervalNode
from random import randint, seed

def normalize(a):
    b = a
    s = sum(a)

    if s == 0:
        for i in a:
            i = 1 / len(a)
    else:
        for i in a:
            i = i / (s + 0.0001)
    return a





def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))




def is_junction(read):
    
    for c in read.cigartuples:
        if c[0]==3:
            return True
    return False





def whichFeature(read,chr):
    find_list_cds=find(read.reference_start, read.reference_end , tree_cds[chr])
    find_list_utr3=find(read.reference_start, read.reference_end , tree_utr3[chr])
    find_list_utr5=find(read.reference_start, read.reference_end , tree_utr5[chr])
    find_list_intron=find(read.reference_start, read.reference_end , tree_geneCoordinates[chr])
    find_list_intergenic=find(read.reference_start, read.reference_end , tree_intergenic[chr])
    

    
    threshold=len(read.query_sequence)*0.75
    
    
    tag_cds=0
    tag_utr3=0
    tag_utr5=0
    max_cds=0
    max_utr3=0
    max_utr5=0
    
    cds_overlap=[]
    utr3_overlap=[]
    utr5_overlap=[]
    
    for i in find_list_cds:
        overlap=getOverlap((read.reference_start,read.reference_end),i)
        if overlap>threshold:
            cds_overlap.append(overlap)
    
    
    for i in find_list_utr3:
        overlap=getOverlap((read.reference_start,read.reference_end),i)
        if overlap>threshold:
            utr3_overlap.append(overlap)
            
    for i in find_list_utr5:
        overlap=getOverlap((read.reference_start,read.reference_end),i)
        if overlap>threshold:
            utr5_overlap.append(overlap)
    

 
    
    if len(cds_overlap)>0:
        tag_cds=1
        max_cds=max(cds_overlap)
    if len(utr3_overlap)>0:
        tag_utr3=1
        max_cds=max(utr3_overlap)
    if len(utr5_overlap)>0:
        tag_utr5=1
        max_cds=max(utr5_overlap)

    
    
    if tag_cds>1 and tag_utr3+tag_utr5>1:
        #print "-------"
        #print "MIXED"
        #print tag_cds,tag_utr3,tag_utr5
        #print max_cds,tag_utr3,tag_utr5
        #print find_list_cds
        #print find_list_utr3
        #print find_list_utr5
        #print read
        #print "-------"
        return 'MIXED'
    elif tag_utr3+tag_utr5>1:
        x=find_list_utr3[0][0]
        y=find_list_utr3[0][1]
        return ('UTR_',geneUTR3[(chr,x,y)])
    elif tag_cds==1:
        
        x=find_list_cds[0][0]
        y=find_list_cds[0][1]
        return ('CDS',geneCDS[(chr,x,y)])
    elif tag_utr3==1:
        x=find_list_utr3[0][0]
        y=find_list_utr3[0][1]
        return ('UTR3',geneUTR3[(chr,x,y)])
    elif tag_utr5==1:
        x=find_list_utr5[0][0]
        y=find_list_utr5[0][1]
        return ('UTR5',geneUTR5[(chr,x,y)])
    elif tag_cds+tag_utr3+tag_utr5==0:
        if len(find_list_intron)>0:
            x=find_list_intron[0][0]
            y=find_list_intron[0][1]
            
            
            return ('INTRON',dictGeneNames[x,y])
        elif len(find_list_intergenic)>0:
            return ('INTERGENIC',("NA","NA"))
        else:
            return ('DEEP',("NA","NA"))

        
        

#------

def find(start, end, tree):
    "Returns a list with the overlapping intervals"
    out = []
    tree.intersect( start, end, lambda x: out.append(x) )
    return [ (x.start, x.end) for x in out ]

'''
    tree = IntervalNode( 5, 20 )
    
    
    
    overlap = find(27, 28 , tree)
    if overlap==[]:
    print "----"
    
    '''


#==============================================================================================================================================
#==============================================================================================================================================
#==============================================================================================================================================
#==============================================================================================================================================
#==============================================================================================================================================
#==============================================================================================================================================


ap = argparse.ArgumentParser()
ap.add_argument('bam', help='sorted bam file with mapped reads')
ap.add_argument('out', help='file to save the number of reads per genome category')
ap.add_argument("--mouse", help="Use mouse genome annotations (NCBIM37). Default is human",action="store_true")
ap.add_argument("--multi", help="Categories all copies of multi-mapped reads, they will have a special flag to futher assign according to transcript abundance ",action="store_true")
args = ap.parse_args()








##About Interval trees  https://www.biostars.org/p/99/

print ("Open bam file",args.bam)
bamfile = pysam.Samfile(args.bam, "rb")



outDir=os.path.dirname(args.out)
if outDir=="":
    outDir=os.getcwd()

chr_list=[]


                
                
if not args.mouse:
    print ("Human annotations are used ...")
    for i in range(1,23):
        chr_list.append(str(i))
    chr_list.append('X')
    chr_list.append('Y')
    chr_list.append('MT')

    utr3_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/UTR3_GRCh37_prepared.bed'
    utr5_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/UTR5_GRCh37_prepared.bed'
    cds_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/CDS_GRCh37_prepared.bed'
    geneCoordinates_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/geneCoordinates_GRCh37.bed'

elif args.mouse:
    print ("Mouse annotations are used ...")
    for i in range(1,20):
        chr_list.append(str(i))
    chr_list.append('X')
    chr_list.append('Y')
    chr_list.append('MT')
    utr3_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/UTR3_NCBIM37_prepared.bed'
    utr5_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/UTR5_NCBIM37_prepared.bed'
    cds_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/CDS_NCBIM37_prepared.bed'
    geneCoordinates_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/geneCoordinates_NCBIM37.bed'



dirOutPerCategory=""


base=os.path.basename(args.bam)
prefix=os.path.splitext(base)[0]






    
    
    
dirOutPerCategory=args.out+"/"
if not os.path.exists(dirOutPerCategory):
        os.makedirs(dirOutPerCategory)
    
print ("Directory to save the results  ", dirOutPerCategory)
    
f_file=dirOutPerCategory+prefix+".genomicFeature"
outFile = open(f_file, 'w' )
outFile.write('readName,chr,category, geneID, geneName, flag_multiMapped')
outFile.write('\n')














#DATA STRUCTURE - per chr
tree_utr3={}
tree_utr5={}
tree_cds={}
tree_geneCoordinates={}
tree_rRNA={}
tree_intergenic={} # +10,000


for chr in chr_list:
    tree_utr3[chr]=IntervalNode(0,0)
    tree_utr5[chr]=IntervalNode(0,0)
    tree_cds[chr]=IntervalNode(0,0)
    tree_geneCoordinates[chr]=IntervalNode(0,0)
    tree_rRNA[chr]=IntervalNode(0,0)
    tree_intergenic[chr]=IntervalNode(0,0)










print ("Load gene annotations ...")


geneUTR3={}

#UTR3
print ("Load",utr3_file)
with open(utr3_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[1])
            y=int(line[2])
            geneID=line[3]
            geneName=line[4]
            geneUTR3[(chr,x,y)]=(geneID,geneName)
            tree_utr3[chr]=tree_utr3[chr].insert( x, y )








#find_list=find(67208778, 67210057 , tree_utr3[chr])

geneUTR5={}

#UTR5
print ("Load",utr5_file)
with open(utr5_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[1])
            y=int(line[2])
            geneID=line[3]
            geneName=line[4]
            geneUTR5[(chr,x,y)]=(geneID,geneName)
            tree_utr5[chr]=tree_utr5[chr].insert( x, y )


#CDS
geneCDS={}


print ("Load",cds_file)
with open(cds_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[1])
            y=int(line[2])
            geneID=line[3]
            geneName=line[4]
            geneCDS[(chr,x,y)]=(geneID,geneName)
            tree_cds[chr]=tree_cds[chr].insert( x, y )






#gene coordinates
nGenes_non_rRNA=0
nGenes_rRNA=0


#['1', 'non-rRNA', 'ENSG00000000457', 'SCYL3', '169818772', '169863408']


print ("Load",geneCoordinates_file)

geneNameSet={}

for c in chr_list:
    geneNameSet[c]=set()

#1,non-rRNA,ENSMUSG00000000544,Gpa33,168060369,168096

dictGeneNames={}

with open(geneCoordinates_file,'r') as f:
    
    reader=csv.reader(f)
    for line in reader:
        chr=line[0]
        if chr in chr_list:
            x=int(line[4])
            y=int(line[5])
            geneID=line[2]
            geneName=line[3]
            geneNameSet[chr].add(geneName)
            
            if line[1]=='non-rRNA':
                nGenes_non_rRNA+=1
                tree_geneCoordinates[chr]=tree_geneCoordinates[chr].insert( x, y )
                dictGeneNames[x,y]=(geneID,geneName)
                
                
                x_10K=x-10000
                y_10K=y+10000
                if x_10K<0:
                    x_10K=0
                tree_intergenic[chr]=tree_intergenic[chr].insert( x_10K, y_10K )
            elif line[1]=='rRNA':
                tree_rRNA[chr]=tree_rRNA[chr].insert( x, y )
                nGenes_rRNA+=1



for c in chr_list:
    print ("Number of genes in the annotations for chr %s = %i" %(c,len(geneNameSet[c])))





#
#======================================================================
#BAM




# the same read can be assigned to multiple genes from different chr


for chr in chr_list:


    #print "Process chr",chr
    for read in bamfile.fetch(chr):
        readName=read.qname



        

        
        flagMulti=0
        if read.get_tag("NH") != 1:
            flagMulti=1
        




        
        if is_junction(read):
            feature=whichFeature(read,chr)
            if flagMulti==0:
                    #1943766_h_0_GCGGATC_GATC,13,CDS,ENSMUSG00000021782,Dlg5,0

                    outFile.write( readName+','+chr + ',' + 'junction' + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )
        
            elif args.multi:
                    outFile.write( readName+','+chr + ',' + 'junction' + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )

                
    
        else:
            feature=whichFeature(read,chr)
            

            if flagMulti==0:
                outFile.write( readName+','+chr + ',' + feature[0] + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )
            elif args.multi:
                outFile.write(readName+','+chr + ',' + feature[0] + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )

outFile.close()

print ("Results are saved to ",f_file)
print ("Done!")



#-===============
# SECOND PART


geneNameSet = {}
geneIDSet = set()

for c in chr_list:
    geneNameSet[c] = set()

# 1,non-rRNA,ENSMUSG00000000544,Gpa33,168060369,168096

dictGeneIds = {}
dict_id2name = {}
dict_id2chr = {}

with open(geneCoordinates_file, 'r') as f:
    reader = csv.reader(f)
    for line in reader:
        chr = line[0]
        if chr in chr_list:
            x = int(line[4])
            y = int(line[5])
            geneID = line[2]
            geneName = line[3]
            geneIDSet.add(geneID)
            dict_id2name[geneID] = geneName
            dict_id2chr[geneID] = chr

print("Number of genes in the annotations = %i" % (len(geneIDSet)))



abundanceGene = {}
mReadsSet = set()
readDict = {}

# 5463462_h_0_TTAAACT_TAGC,1,CDS,ENSMUSG00000051951,Xkr4,0


# updated 01/14/2017
# .genomicFeature files are parsed and we ramdonly select were multi-mapped reads is mapped
# We select ramdomly:
# 1. read is mapped to multiple location of the same gene
# 2. read is mapped otside of the gene . Or mapped into a gene and outsde of the gene
# We use relative frequency of the gene (calculated based on uniqelly mapped reads)
# 1. read is mapped to multiple location, all those locations are within different genes






for g in geneIDSet:
    abundanceGene[g] = 0

# f_file with reads mapped to multiple genes, i.e. multiple lines per gene





# f_file2 with reads, where one location was choosen for read randomly
f_file2 = dirOutPerCategory + prefix + ".genomicFeature2"

out = open(f_file2, 'w')
out.write("readName,chr,category,geneID,geneName,flag_multiMapped")
out.write("\n")

with open(f_file, 'r') as f:
    reader = csv.reader(f)
    next(reader, None)

    for line in reader:
        read = line[0]
        geneID = line[3]
        flagM = int(line[5])
        readName = line[0]
        mReadsSet.add(readName)

for r in mReadsSet:
    readDict[r] = []

reads_unique_not_used_counts = 0
# print "Reading",f_file
with open(f_file, 'r') as f:
    reader = csv.reader(f)
    next(reader, None)
    for line in reader:
        category = line[2]
        geneID = line[3]
        flagM = int(line[5])
        readName = line[0]
        readDict[readName].append(line)
        # abundanceGene - we need this to calculate
        if geneID != "NA" and flagM == 0 and (category == "CDS" or category=="junction"):
            abundanceGene[geneID] += 1
        elif flagM == 0:
            reads_unique_not_used_counts += 1

print("sum(abundanceGene.values()),reads_unique_not_used_counts", sum(abundanceGene.values()),
      reads_unique_not_used_counts)



count_intergenic=0

for r in mReadsSet:

    # read mapped to multiple locations of the same gene

    weights = []
    weights[:] = []

    if len(readDict[r]) == 1:  # some reads were multi-mapped, but than got collapsed and now have flag 1 but are present in a single copy in [genomicFeature] file
        out.write(readDict[r][0][0] + "," + readDict[r][0][1] + "," + readDict[r][0][2] + "," + readDict[r][0][3] + "," +readDict[r][0][4] + "," + readDict[r][0][5])
        out.write("\n")


    if len(readDict[r]) > 1:  #read is multi-mapped




        # print "==============="
        # print readDict[r], "MULTI-MAPPED"

        genes = set()
        genes.clear()

        readsList=[]
        readsList[:]=[]


        for g in readDict[r]:
            if g[2]=="CDS" or g[2]=="junction":
                if abundanceGene[g[3]]!=0:
                    readsList.append(g)




        irand = 0
        if len(readsList)==1: #only one read from CDS left
            out.write(
                readsList[irand][0] + "," + readsList[irand][1] + "," + readsList[irand][2] + "," + readsList[irand][
                    3] + "," + readsList[irand][4] + "," + readsList[irand][5])
            out.write("\n")
        elif len(readsList)>1:
            listIndex = list(range(len(readsList)))
            for g in readsList:
                weights.append(abundanceGene[g[3]])
            print (weights)
            norm = [float(i) / sum(weights) for i in weights]
            weights = norm
            irand = choice(listIndex, p=weights)
            print ("irand",irand)
            out.write(readsList[irand][0] + "," + readsList[irand][1] + "," + readsList[irand][2] + "," +readsList[irand][3] + "," + readsList[irand][4] + "," + readsList[irand][5])
            out.write("\n")












        # make sure to clean everething, as we fo to new chromosome

out.close()

f_summary = dirOutPerCategory + prefix + "_summary_per_feature.csv"
fileS = open(f_summary, "w")
fileS.write("nJunction,nCDS,nUTR3,nUTR5,nUTR_,nIntron")
fileS.write("\n")

geneIDSetCurrent = set()
dict = {}
for g in geneIDSet:
    dict[g] = [0, 0, 0, 0, 0, 0]



InterGenic_set=set()

f_fileOut = dirOutPerCategory + prefix + ".perGeneSummary"
with open(f_file2, 'r') as f:
    reader = csv.reader(f)
    for line in reader:

        if len(line) == 6:  # fix it later
            read = line[0]
            geneID = line[3]

            geneIDSetCurrent.add(geneID)
            category = line[2]

            if category == "junction" and geneID != "NA":
                dict[geneID][0] += 1
            elif (category == "CDS"):
                dict[geneID][1] += 1
            elif (category == "UTR3"):
                dict[geneID][2] += 1
            elif (category == "UTR5"):
                dict[geneID][3] += 1
            elif category == "UTR_":
                dict[geneID][4] += 1
            elif category == "INTRON":
                dict[geneID][5] += 1

        else:
            print("Warning:", line)

outfile = open(f_fileOut, 'w')
outfile.write("ID,geneName,chr,nJunction,nCDS,nUTR3,nUTR5,nUTR_,nIntron,gene.count,gene.count.unique \n")

nJunction = 0
nCDS = 0
nUTR3 = 0
nUTR5 = 0
nUTR_ = 0
nIntron = 0

for g in geneIDSet:
    if g in geneIDSetCurrent:
        nJunction += dict[g][0]
        nCDS += dict[g][1]
        nUTR3 += dict[g][2]
        nUTR5 += dict[g][3]
        nUTR_ += dict[g][4]
        nIntron += dict[g][5]
        outfile.write(g + "," + dict_id2name[g] + "," + dict_id2chr[g] + "," + str(dict[g][0]) + "," + str(
            dict[g][1]) + "," + str(dict[g][2]) + "," + str(dict[g][3]) + "," + str(dict[g][4]) + "," + str(
            dict[g][5]) + "," + str(sum(dict[g])) + "," + str(abundanceGene[g]) + "\n")
    else:
        outfile.write(g + "," + dict_id2name[g] + "," + dict_id2chr[g] + "," + "0,0,0,0,0,0\n")

outfile.close()

fileS.write(str(nJunction) + "," + str(nCDS) + "," + str(nUTR3) + "," + str(nUTR5) + "," + str(nUTR_) + "," + str(
    nIntron))
fileS.write("\n")
fileS.close()

print("Summary per category is here", f_summary)
print("Summary per gene is here", f_fileOut)
print("Done!")




