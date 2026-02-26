#!/bin/bash
SECONDS=0

# FASTQC
run_fastqc (){
input_dir="$1"
output_dir="$2"
filename_pattern="$3"
conda_env="Fast_QC"
echo "Ejecutando FastQC en archivos de $input_dir con patrón $filename_pattern..."
mkdir -p "$output_dir"
conda activate "$conda_env"
find "$input_dir" -type f -name "$filename_pattern" | while read -r file; do
    echo "Procesando $file..."
    fastqc "$file" -o "$output_dir"
done
echo "FastQC finalizado. Resultados en $output_dir"
}


#MULTIQC4FASTQC
run_multiqc (){
  input_dir="$1"
  output_dir="$2"
  filename_pattern="$3"
# Nombre del entorno conda donde está MultiQC (ajusta si usas otro nombre)
  conda_env="Multi_QC"

  echo "Ejecutando MultiQC en archivos de $input_dir con patrón $filename_pattern..."
  mkdir -p "$output_dir"

  # Activar conda (asume que 'conda' es accesible y configurado en el sistema)
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$conda_env"

  # Ejecutar MultiQC sobre los archivos que coincidan con el patrón
  multiqc "$input_dir"/$filename_pattern -o "$output_dir"

  # Desactivar el entorno conda
  conda deactivate

  echo "MultiQC finalizado. Reporte en $output_dir"
}

 #TRIMMOMATIC4CROPPINGADAPTERS
#   $3 = R1_pattern (por defecto "*_R1.001.fastq")
#   $4 = threads (por defecto 8)
#   $6 = adapters_path (por defecto "/path/to/TruSeq3-PE-2.fa") -- ajústalo a la ruta real en tu sistema

# se corre./RunPipeline.sh run_trimmomatic \ "/Users/JK/Desktop/JKViromedata/Aa" \ "/Users/JK/D
# se corre no valores default./RunPipeline.sh run_trimmomatic "/Users/JK/Desktop/JKViromedata/Aa" "/Users/JK/Desktop/SecRef4TrimmAdapters" "*_R1.001.fastq" \10 \trimmomatic_env "/Users/Parsimony/miniconda3/envs/trimmomatic/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa"
run_trimmomatic (){
  input_dir="$1"
  output_dir="$2"
  r1_pattern="$3"
  threads="${4:-8}"
  conda_env="${5:-trimmomatic_env}"
  echo "Ejecutando Trimmomatic en $input_dir con patrón $r1_pattern..."
  mkdir -p "$output_dir"
  # Activar conda y entorno
  export _JAVA_OPTIONS="-Xmx8G"
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$conda_env"
  # Loop simple sobre archivos R1 que coincidan con el patrón
for R1 in "$input_dir" ; do
    R2="${R1/_R1/_R2}"

    sample_name=$(basename "$R1" _R1.fastq)
    echo "Sample: $sample_name"
    echo "R1: $R1"
    echo "R2: $R2"
    echo "--------------"

    echo "Procesando pareja: ${BASENAME} (R1: $R1, R2: $R2)"
    trimmomatic PE -threads "$threads" -phred33 \
      "$R1" "$R2" \
      "$output_dir/${BASENAME}_R1_paired.fastq.gz" \
      "$output_dir/${BASENAME}_R1_unpaired.fastq.gz" \
      "$output_dir/${BASENAME}_R2_paired.fastq.gz" \
      "$output_dir/${BASENAME}_R2_unpaired.fastq.gz" \
      ILLUMINACLIP:/opt/anaconda3/envs/trimmomatic_env/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10:2:keepBothReads \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:35
  done
  conda deactivate

  echo "Trimmomatic finalizado. Salida en $output_dir"
}
 
# RUNSTAR4INDEXING
 # Posicionales:
#   $5 = conda_env          (opcional, por defecto "star_env")
#   $6 = sjdbOverhang       (opcional, por defecto 109)
#   $7 = genomeSAindexNbases (opcional, por defecto 14) genomas pequeños
# ./RunPipeline.sh run_star_index "/Users/JK/Desktop/Aae_refgenome/GCF_002204515.2_AaegL5.0_genomic.fna" "/Users/JK/Desktop/Aae_refgenome/genomic.gtf" "/Users/JK/Desktop/STARgenomeIndexAae" \4 \star_env \109 \14
run_star_index (){
  genome_fasta="$1"
  gtf_file="$2"
  genome_dir="$3"
  threads="${4:-8}"
  conda_env="${5:-star_env}"

  echo "Ejecutando STAR genomeGenerate..."
  mkdir -p "$genome_dir"
  chmod -R 777 "$genome_dir"
  chmod 644 "$gtf_file"

  # Activar conda (se asume que conda está disponible)
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$conda_env"
  STAR --runThreadN "$threads" \
       --runMode genomeGenerate \
       --genomeDir "$genome_dir" \
       --genomeFastaFiles "$genome_fasta" \
       --sjdbGTFfile "$gtf_file" \

  conda deactivate

  echo "STAR genome indexing finalizado. Output: $genome_dir"
}

# RUNSTAR4MAPPING

# Posicionales:
#   $1 = input_dir        (directorio con los *_R1_paired.fastq.gz)
#   $2 = genome_dir       (directorio del índice STAR)
#   $3 = gtf_file         (ruta al .gtf)
#   $4 = outdir           (directorio de salida para los alineamientos)
#   $5 = r1_pattern       (patrón para R1, por defecto "*R1.001.fastq")
#   $6 = threads          (por defecto 8)
#   $7 = conda_env        (por defecto "mapping")
run_star_align (){
  input_dir="$1"
  genome_dir="$2"
  gtf_file="$3"
  outdir="$4"
  r1_pattern="$5"
  threads="${6:-8}"
  conda_env="${7:-mapping}"

  echo "Ejecutando STAR align en $input_dir con patrón $r1_pattern..."

  mkdir -p "$outdir"

  # Activar conda (se asume que 'conda' está disponible)
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "$conda_env"

 for R1 in "$input_dir"/*"$r1_pattern"
  do
    # Chequeo de existencia
    if [[ ! -f "$R1" ]]; then
      echo "Archivo $R1 no encontrado, saltando..."
      continue
    fi
    # Quita el sufijo _R1.fastq (ajusta según tu nomenclatura real)
    BASENAME=$(basename "$R1" "$r1_pattern")
    # Genera el nombre de R2 cambiando "_R1.fastq" por "_R2.fastq"
    R2="${input_dir}/${BASENAME}_R2.fastq"

    # Chequea si R2 realmente existe antes de continuar
    if [[ ! -f "$R2" ]]; then
      echo "Pareja no encontrada para $R1, falta $R2. Saltando..."
      continue
    fi
    sample_names=$(basename "$R1" "R1.fastq")
    echo "Procesando muestra: $sample_names"
    STAR --runThreadN "$threads" \
         --runMode alignReads \
         --genomeDir "$genome_dir" \
         --readFilesIn "$R1" "$R2" \
         --outFileNamePrefix "${outdir}${sample_names}_"
  done

  conda deactivate

  echo "STAR align finalizado. Salida en: $outdir"
}



# se corre asi: ./RunPipeline.sh run_fastqc "/Users/JK/Desktop/JKViromedata/Aa" "/Users/JK/Desktop/JKViromedata" "*.fastq"
# ./RunPipeline.sh run_star_index "/Users/JK/Desktop/Aae_refgenome/GCF_002204515.2_AaegL5.0_genomic.fna" "/Users/JK/Desktop/Aae_refgenome/genomic.gtf" "/Users/JK/Desktop/STARgenomeIndexAae" \4 \star_env \109 \14. Ejemplo practico
# Selecciona la función según el primer argumento
# Dispatcher sencillo para invocar las funciones
case "$1" in
  run_fastqc)
    run_fastqc "$2" "$3" "$4"
    ;;
  run_multiqc)
    run_multiqc "$2" "$3" "$4"
    ;;
    run_trimmomatic)
    run_trimmomatic "$2" "$3" "$4" "$5" "$6"
    ;;
    run_star_index)
    run_star_index "$2" "$3" "$4" "$5" "$6" "$7" "$8"
    ;;
    run_star_align)
    run_star_align "$2" "$3" "$4" "$5" "$6" "$7" "$8"
    ;;
  *)
  echo "Uso: $0 {run_fastqc|run_multiqc|run_trimmomatic|run_star_index|run_star_align} <input_dir> <output_dir> \"<filename_pattern>\"" #{run_lafuncionquequieras}
esac

    duration=$SECONDS
    echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."