import pysam
import sys
import csv
import os
import argparse


def normalize(a):
    b=a
    s=sum(a)
    
    if s==0:
        for i in a:
            i=1/len(a)
    else:
        for i in a:
            i=i/(s+0.0001)
    return a


ap = argparse.ArgumentParser()
ap.add_argument('dir', help='perCategory directory generated by running gprofile.py')
ap.add_argument("--mouse", help="Use mouse genome annotations (NCBIM37). Default is human",action="store_true")

args = ap.parse_args()





dirOutPerCategory=args.dir





chr_list=[]


                
                
if not args.mouse:
    print "Human annotations are used ..."
    for i in range(1,23):
        chr_list.append(str(i))
    chr_list.append('X')
    chr_list.append('Y')
    chr_list.append('MT')

    geneCoordinates_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/human/bedPrepared/geneCoordinates_GRCh37.bed'

elif args.mouse:
    print "Mouse annotations are used ..."
    for i in range(1,20):
        chr_list.append(str(i))
    chr_list.append('X')
    chr_list.append('Y')
    chr_list.append('MT')
    geneCoordinates_file=os.path.dirname(os.path.realpath(__file__))+'/annotations/mouse/bedPrepared/geneCoordinates_NCBIM37.bed'











#gene coordinates
nGenes_non_rRNA=0
nGenes_rRNA=0


#['1', 'non-rRNA', 'ENSG00000000457', 'SCYL3', '169818772', '169863408']


print "Load",geneCoordinates_file

geneNameSet={}
geneIDSet={}

for c in chr_list:
    geneNameSet[c]=set()
    geneIDSet[c]=set()

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
            geneIDSet[chr].add(geneID)
            dictGeneNames[geneName]=geneID


for c in chr_list:
    print "Number of genes in the annotations for chr %s = %i" %(c,len(geneNameSet[c]))

print "Look into ", args.dir




for f in os.listdir(dirOutPerCategory):
    if f.endswith(".genomicFeature"):
        base=os.path.basename(f)
        if "MT.genomicFeature" in f:
           prefix=os.path.splitext(base)[0].split('.MT.genomicFeature')[0].replace('.MT','') # not sure why we need replace. But somehow split works wired 09/01/2017





abundanceGene={}
mReadsSet=set()
readDict={}





#5463462_h_0_TTAAACT_TAGC,1,CDS,ENSMUSG00000051951,Xkr4,0


#updated 01/14/2017
# .genomicFeature files are parsed and we ramdonly select were multi-mapped reads is mapped
# We select ramdomly:
# 1. read is mapped to multiple location of the same gene
# 2. read is mapped otside of the gene . Or mapped into a gene and outsde of the gene
# We use relative frequency of the gene (calculated based on uniqelly mapped reads)
# 1. read is mapped to multiple location, all those locations are within different genes


reads_set_MT=set()

for chr in chr_list:
    
    for g in geneIDSet[chr]:
        abundanceGene[g]=0
    
    #f_file with reads mapped to multiple genes, i.e. multiple lines per gene
    
    print "dirOutPerCategory",dirOutPerCategory
    print "prefix",prefix
    print "chr",chr
    
    f_file=dirOutPerCategory+prefix+"."+chr+".genomicFeature"

    #f_file2 with reads, where one location was choosen for read randomly
    f_file2=dirOutPerCategory+prefix+"."+chr+".genomicFeature2"


    out = open(f_file2, 'w')



    with open(f_file,'r') as f:
        reader=csv.reader(f)
        
        reader.next()
        
        for line in reader:
            read=line[0]
            geneID=line[3]
            flagM=int(line[5])
            readName=line[0]
            mReadsSet.add(readName)

            if chr=="MT":
                reads_set_MT.add(read)


    for r in mReadsSet:
        readDict[r]=[]



    print "Reading",f_file
    with open(f_file,'r') as f:
        reader=csv.reader(f)
        
        reader.next()
        
        for line in reader:
            
            geneID=line[3]
            flagM=int(line[5])
            readName=line[0]
            
           
            
        
        
        
            readDict[readName].append(line)
                        
            if flagM==1 and readName not in mReadsSet:
                abundanceGene[geneID]=0



            if geneID!="NA" and flagM==0:
                abundanceGene[geneID]+=1



    



    from numpy.random import choice
    import random


    for r in mReadsSet:
        
      
        
        #read mapped to multiple locations of the same gene
        
        
        weights = []
        weights[:]=[]
        
        
        if len(readDict[r])==1: # some reads were multi-mapped, but than got collapsed and now have flag 1 but are present in a single copy in [genomicFeature] file
            out.write(readDict[r][0][0]+","+readDict[r][0][1]+","+readDict[r][0][2]+","+readDict[r][0][3]+","+readDict[r][0][4]+","+readDict[r][0][5])
            out.write("\n")
        
        if len(readDict[r])>1:#read is multi-mapped
            
            
            
            #print "==============="
            #print readDict[r], "MULTI-MAPPED"
            
            genes=set()
            genes.clear()
            for g in readDict[r]:
                genes.add(g[3])
            
            irand=0
            
            if len(genes)==1 or "NA" in genes: #read mapped to different location within the same gene OR genes vs intergenic
                
                #print "-->(1) Is mapped to different location within the same gene"
                
                
                irand = random.randrange(0, len(readDict[r]))
                
                
                
            
            else:
                
                
                listIndex=list(range(len(readDict[r])))
                
                #print "-->(2)Is mapped to different genes"

                
                for g in readDict[r]:
                    weights.append(abundanceGene[g[3]])
                
                
                if sum(weights)==0.0: # in case read belongs to genrs with no uniq reads. As the results the freq of those genes is 0. And we can not do assigment based on freq of genes
                    #print "Sum is 0!"
                    
                    norm=[float(1.0/len(weights))] * len(weights)
                else:
                    norm = [float(i)/sum(weights) for i in weights]

                weights=norm
                
                
                irand=choice(listIndex, p=weights)
            
            #print "We choose",irand,readDict[r][irand]
            
            
            
            out.write(readDict[r][irand][0]+","+readDict[r][irand][1]+","+readDict[r][irand][2]+","+readDict[r][irand][3]+","+readDict[r][irand][4]+","+readDict[r][irand][5])
            out.write("\n")


            #make sure to clean everething, as we fo to new chromosome






f_summary=dirOutPerCategory+prefix+"_summary_per_feature.csv"
fileS=open(f_summary,"w")
fileS.write("chr,nJunction,nCDS,nUTR3,nUTR5,nUTR_,nIntron")
fileS.write("\n")

nJunction_g=0
nCDS_g=0
nUTR3_g=0
nUTR5_g=0
nUTR_g=0
nIntron_g=0



#new stuff added by Sarah request
#updated before releasing UMI-Reducer 01/11/2017

intergenicReads=set()

for chr in chr_list:
    
    
    
    
    
    if not os.path.exists(dirOutPerCategory+"/perGeneSummary/"):
        os.makedirs(dirOutPerCategory+"/perGeneSummary/")
    
    
    
    f_file=dirOutPerCategory+prefix+"."+chr+".genomicFeature2"
    f_fileOut=dirOutPerCategory+"/perGeneSummary/"+prefix+"."+chr+".perGeneSummary"
    
    print "Processing ",f_file





    genes=set()

    with open(f_file,'r') as f:
        
        reader=csv.reader(f)
        for line in reader:
            
            if len(line)==6:
                genes.add(line[4])


    #nJunction,nCDS,nUTR3,nUTR5,nUTR_,nIntron

    dict={}
    for g in genes:
        dict[g]=[0,0,0,0,0,0]

    geneNameSetCurrent=set()
    geneNameSetCurrent.clear()

    with open(f_file,'r') as f:
        
        reader=csv.reader(f)
        for line in reader:
            
            if len(line)==6:
                read=line[0]
                geneName=line[4]
                geneNameSetCurrent.add(geneName)
                
                category=line[2]
                
               
                
                
                if (category =="junction"):
                    dict[geneName][0]+=1
                elif (category =="CDS"):
                    dict[geneName][1]+=1
                elif (category =="UTR3"):
                    dict[geneName][2]+=1
                elif (category =="UTR5"):
                    dict[geneName][3]+=1
                elif category=="UTR_":
                    dict[geneName][4]+=1
                elif category=="INTRON":
                    dict[geneName][5]+=1
                else:
                    intergenicReads.add(read)

    outfile = open(f_fileOut, 'w' )


    outfile.write("geneName,chr,nJunction,nCDS,nUTR3,nUTR5,nUTR_,nIntron\n")

    nJunction=0
    nCDS=0
    nUTR3=0
    nUTR5=0
    nUTR_=0
    nIntron=0
    

    


    for g in geneNameSet[chr]:
        if g in geneNameSetCurrent:
            value=dict[g]
            
            geneId=dictGeneNames[g]
            
            
            
            nJunction+=value[0]
            nCDS+=value[1]
            nUTR3+=value[2]
            nUTR5+=value[3]
            nUTR_+=value[4]
            nIntron+=value[5]
            
            
            
            outfile.write(g+","+chr+","+str(value[0])+","+str(value[1])+","+str(value[2])+","+str(value[3])+","+str(value[4])+","+str(value[5])+","+str(sum(value))+","+str(abundanceGene[geneId])+"\n")
        else:
            outfile.write(g+","+chr+","+"0,0,0,0,0,0\n")

    outfile.close()


    nJunction_g+=nJunction
    nCDS_g+=nCDS
    nUTR3_g+=nUTR3
    nUTR5_g+=nUTR5
    nUTR_g+=nUTR_
    nIntron_g+=nIntron


    if chr!="MT":
        fileS.write(chr+","+str(nJunction)+","+str(nCDS)+","+str(nUTR3)+","+str(nUTR5)+","+str(nUTR_)+","+str(nIntron))
        fileS.write("\n")


fileS.write("Total,"+str(nJunction_g)+","+str(nCDS_g)+","+str(nUTR3_g)+","+str(nUTR5_g)+","+str(nUTR_g)+","+str(nIntron_g))
fileS.write("\n")








