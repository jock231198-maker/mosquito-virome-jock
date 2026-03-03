sudo cat M11_S3_L001_R1_001.fastq.gz M11_S3_L002_R1_001.fastq.gz > M11_S3_R1_001.fastq.gz
cat M11_S3_L001_R1_001.fastq.gz M11_S3_L002_R1_001.fastq.gz > /Users/JK/Desktop/M11_S3_R1_001.fastq.gz
cat M11_S3_L001_R2_001.fastq.gz M11_S3_L002_R2_001.fastq.gz > /Users/JK/Desktop/M11_S3_R2_001.fastq.gz

STAR \
--runMode alignReads \
--runThreadN 16 \
--genomeDir /Users/JK/Desktop/STARresults/Star_Index \
--readFilesIn /Users/JK/Desktop/M11_S3_R1_001.paired.fastq /Users/JK/Desktop/M11_S3_R2_001.paired.fastq \
--outFileNamePrefix /Users/JK/Desktop/STARresults/Star_Map/M1_reads \

STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeDir /Users/JK/Desktop/STARresults/Star_index_albo \
  --genomeFastaFiles /Users/JK/Desktop/GCF_035046485.1_AalbF5_genomic.fna \
  --sjdbGTFfile /Users/JK/Desktop/genomic.gtf \
  --sjdbOverhang 99 \

  STAR \
--runMode alignReads \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 8 \
--outReadsUnmapped Fastx \
--outMultimapperOrder Random \
--outWigType wiggle \
--genomeDir /Users/JK/Desktop/STARresults/Star_index_albo \
--readFilesIn /Users/JK/Desktop/M11_S3_R1_001.fastq /Users/JK/Desktop/M11_S3_R2_001.fastq \
--outFileNamePrefix /Users/JK/Desktop/STARresults/Star_Map/M1_reads \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--outFilterMultimapNmax 10 
