#!/bin/bash
# STAR align
{                                 
     OUTDIR="/Users/JK/Desktop/STARresults/STARalign/StarSimple/"
     mkdir -p "$OUTDIR"
     for R1 in /Users/JK/Desktop/JKViromedata/M1_M11/*R1.fastq ; do
    R2="${R1/_R1/_R2}"
    sample_name=$(basename "$R1" _R1.fastq.gz)

    echo "Sample: $sample_name"
    echo "R1: $R1"
    echo "R2: $R2"
  STAR \
        --runThreadN 4 \
        --runMode alignReads \
        --genomeDir /Users/JK/Desktop/STARresults/StarIndexAae \
        --readFilesIn "$R1" "$R2" \
        --outFileNamePrefix "${OUTDIR}${sample_names}_" \
      
done
}
conda deactivate

