# Download all RefSeq viral genomes
datasets download virus genome taxon Viruses \
  --refseq \
  --filename viral_refseq.zip
  
#unzip and concatenate all fasta files into one
unzip viral_refseq.zip
cat ncbi_dataset/data/**/*.fna > viral_refseq_all.fasta