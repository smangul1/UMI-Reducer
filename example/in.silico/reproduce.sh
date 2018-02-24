
/u/home/s/svandrie/code/import/hisat2-2.1.0/hisat2-build myReferenceGenome.fasta myReferenceGenome

/u/home/s/svandrie/code/import/hisat2-2.1.0/hisat2  myReferenceGenome -U reads.simulated.fastq  | samtools view -F 4 -bS - >reads.simulated.bam  

samtools sort reads.simulated.bam >reads.simulated.sort.bam
samtools index reads.simulated.sort.bam

