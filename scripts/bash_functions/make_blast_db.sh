# Build the nucleotide database
makeblastdb \
  -in viral_refseq_all.fasta \
  -dbtype nucl \
  -title "Viral_RefSeq" \
  -out blast_db/viral_refseq_db \
  -parse_seqids       # allows fetching sequences by ID later