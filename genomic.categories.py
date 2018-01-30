import pysam
import sys
import csv
import os
import argparse



from quicksect import IntervalNode
from random import randint, seed


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
ap.add_argument("--perCategory", help="reports the assigment for each read. A separate file per chromosome will be created. In case --perCategory was choosen this will be directory to sabve results",action="store_true")
ap.add_argument("--mouse", help="Use mouse genome annotations (NCBIM37). Default is human",action="store_true")
ap.add_argument("--multi", help="Categories all copies of multi-mapped reads, they will have a special flag to futher assign according to transcript abundance ",action="store_true")
args = ap.parse_args()








##About Interval trees  https://www.biostars.org/p/99/

print ("Open bam file",args.bam)
bamfile = pysam.AlignmentFile(args.bam, "rb")



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

#-------perCategory-------
if args.perCategory:
    
    
    
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
        readName=read.query_name



        

        
        flagMulti=0
        if read.get_tag("NH") != 1:
            flagMulti=1
        
        
        
        if is_junction(read):
            feature=whichFeature(read,chr)
            if flagMulti==0:
                    #1943766_h_0_GCGGATC_GATC,13,CDS,ENSMUSG00000021782,Dlg5,0
                    if args.perCategory:
                        outFile.write( readName+','+chr + ',' + 'junction' + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )
        
            elif args.multi:
                    if args.perCategory:
                        outFile.write( readName+','+chr + ',' + 'junction' + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )

                
    
        else:
            feature=whichFeature(read,chr)
            
            if args.perCategory:
                if flagMulti==0:
                    outFile.write( readName+','+chr + ',' + feature[0] + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )
                elif args.multi:
                    outFile.write( readName+','+chr + ',' + feature[0] + ',' + feature[1][0] + ',' + feature[1][1] + ',' + str(flagMulti)+'\n' )

outFile.close()

print "Results are saved to ",f_file
print "Done!"



