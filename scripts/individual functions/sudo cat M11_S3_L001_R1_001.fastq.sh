$Local runs
sudo cat M11_S3_L001_R1_001.fastq.gz M11_S3_L002_R1_001.fastq.gz > M11_S3_R1_001.fastq.gz
cat M11_S3_L001_R1_001.fastq.gz M11_S3_L002_R1_001.fastq.gz > /Users/JK/Desktop/M11_S3_R1_001.fastq.gz
cat M11_S3_L001_R2_001.fastq.gz M11_S3_L002_R2_001.fastq.gz > /Users/JK/Desktop/M11_S3_R2_001.fastq.gz


STAR \
--runMode alignReads \
--runThreadN 4 \
--genomeDir /Users/JK/Desktop/STARresults/Star_Index \
--readFilesIn /Users/JK/Desktop/M11_S3_R1_001.paired.fastq /Users/JK/Desktop/M11_S3_R2_001.paired.fastq \
--outFileNamePrefix /Users/JK/Desktop/STARresults/Star_Map/M1_reads 

STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir /Users/JK/Desktop/STARresults/Star_index_albo \
  --genomeFastaFiles /Users/JK/Desktop/GCF_035046485.1_AalbF5_genomic.fna \
  --sjdbGTFfile /Users/JK/Desktop/genomic.gtf \

  STAR \
--runMode alignReads \
--outSAMattributes All \
--outSAMtype BAM SortedByCoordinate \
--runThreadN 4 \
--outReadsUnmapped Fastx \
--outMultimapperOrder Random \
--outWigType wiggle \
--genomeDir /Users/JK/Desktop/STARresults/Star_index_albo \
--readFilesIn /Users/JK/Desktop/M11_S3_R1_001.fastq /Users/JK/Desktop/M11_S3_R2_001.fastq \
--outFileNamePrefix /Users/JK/Desktop/STARresults/Star_Map/M1_reads \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--outFilterMultimapNmax 10 

STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir /Users/JK/Desktop/STARresults/Star_Index \
  --genomeFastaFiles /Users/JK/Desktop/AaaGenRef/GCF_002204515.2_AaegL5.0_genomic.fna \
  --sjdbGTFfile /Users/JK/Desktop/AaaGenRef/GCF_002204515.2_AaegL5.0_genomic.gtf 

  STAR \
  --runThreadN 4 \
  --runMode genomeGenerate \
  --genomeDir /Users/JK/Desktop/STARresults/Star_Index \
  --genomeFastaFiles /Users/JK/Desktop/AaaGenRef/GCF_002204515.2_AaegL5.0_genomic.fna 
  
mkdir chmod -m 777 StarAlign
  STAR \
--runMode alignReads \
--runThreadN 4 \
--genomeDir /Users/JK/StarIndex \
--readFilesIn /Users/JK/Desktop/M11_S3_R1_001.paired.fastq /Users/JK/Desktop/M11_S3_R2_001.paired.fastq \
--outFileNamePrefix /Users/JK/StarIndex/StarAlign

make STARforMacStatic CXX=/opt/homebrew/Cellar/gcc/