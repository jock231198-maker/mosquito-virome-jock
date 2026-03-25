bowtie2 
#Indexing
# Building a small index
#comandtobuild /path/to/genome/fna /path/to/index/directory
bowtie2-build /Users/JK/Desktop/GCF_035046485.1_AalbF5_genomic.fna /Users/JK/Desktop/Bow/IndexBow/Aea_index

# Building a large index
bowtie2-build --large-index example/reference/lambda_virus.fa example/index/lambda_virus

#Aligning
# Aligning unpaired reads
bowtie2 -x example/index/lambda_virus -U example/reads/longreads.fq


# Aligning paired reads
#bowtie2 -x(permission to execute) path/to/index/directory -1 path/to/read1 -2 path/to/read2 -S /path/to/outputfolder/namefile.sam
bowtie2 -x /Users/JK/Desktop/Bow/IndexBow -1 /Users/JK/Desktop/M11_S3_R1_001.paired.fastq -2 /Users/JK/Desktop/M11_S3_R2_001.paired.fastq -S /Users/JK/Desktop/Bow/MapBow/M11_S3_mapped.sam --threads 4
