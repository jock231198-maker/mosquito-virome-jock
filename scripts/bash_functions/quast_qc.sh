
mkdir -p /Users/JK/Desktop/quast/M9_S16
quast.py /Users/JK/Desktop/spades/Assembly/trimm_cat_fastq_unmapped/M9_S16/contigs.fasta \
    -o /Users/JK/Desktop/quast/M9_S16 \
    --rna-finding \
    --fast \
    --threads 8 \
    --min-contig 200 