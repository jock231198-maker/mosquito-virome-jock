blastn \
  -query /Users/JK/Desktop/spades/Assembly/trimm_cat_fastq_unmapped/M1_S1/contigs.fasta \
  -db /Users/JK/Desktop/ncbi_database/blast_db/viral_refseq_db\
  -out /Users/JK/Desktop/blast_results/blast_results.txt \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -evalue 1e-5 \
  -num_threads 8 \
  -max_target_seqs 5