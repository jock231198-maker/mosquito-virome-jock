# STAR_index
chmod 644 /Users/JK/Desktop/Aae_refgenome/genomic.gtf #gtf con permisos

STAR \
{
mkdir -p /Users/JK/Desktop/STARgenomeIndexAae 
chmod 755 /Users/JK/Desktop/STARgenomeIndexAae
#crea un directorio de salida con permismos
conda activate star_env
STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir /Users/JK/Desktop/STARresults/Star_Index \
  --genomeFastaFiles /Users/JK/Desktop/AaaGenRef/GCF_002204515.2_AaegL5.0_genomic.fna \
  --sjdbGTFfile /Users/JK/Desktop/AaaGenRef/GCF_002204515.2_AaegL5.0_genomic.gtf 

conda deactivate


# Posicionales:
#   $5 = conda_env          (opcional, por defecto "star_env")
#   $6 = sjdbOverhang       (opcional, por defecto 109)
#   $7 = genomeSAindexNbases (opcional, por defecto 14) genomas pequeños
run_star_index (){
  genome_fasta="$1"
  gtf_file="$2"
  genome_dir="$3"
  threads="${4:-8}"
  conda_env="${5:-star_env}"
  sjdbOverhang="${6:-109}"
  genomeSAindexNbases="${7:-14}"
  echo "Ejecutando STAR genomeGenerate..."
  mkdir -p "$genome_dir"
  chmod 755 "$genome_dir"
  chmod 644 "$gtf_file"

  # Activar conda (se asume que conda está disponible)
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$conda_env"
  STAR --runThreadN "$threads" \
       --runMode genomeGenerate \
       --genomeDir "$genome_dir" \
       --genomeFastaFiles "$genome_fasta" \
       --sjdbGTFfile "$gtf_file" \
       --sjdbGTFtagExonParentTranscript transcript_id \
       --sjdbOverhang "$sjdbOverhang" \
       --genomeSAindexNbases "$genomeSAindexNbases"

  conda deactivate

  echo "STAR genome indexing finalizado. Output: $genome_dir"
}
