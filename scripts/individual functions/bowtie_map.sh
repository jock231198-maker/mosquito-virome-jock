#!/bin/bash
#bowtiemap
SECONDS=0
eval "$(conda shell.bash hook)"
conda activate bow # Inicializa conda # Activar el entorno con Bowtie2 instalado
#runlike this ./bowtie_map.sh "/directory/to/my/viromedata/"
input_dir="$1"
result_from="$2"
index_dir="/Users/JK/Desktop/bowtie2/IndexBow/aedes_aegypti_smallindex/aedes_aegypti_smallindex" # Ruta al índice de Bowtie2
find "$input_dir" -type f -name "*_R1_001_paired.fastq.gz" | while read -r file; do
sample=$(basename "$file" | cut -d'_' -f1-2)
outdir="/Users/JK/Desktop/bowtie2/MapBow/${result_from}/"
mkdir -p "$outdir"
echo "mapping $sample"

read1="$input_dir/${sample}_R1_001_paired.fastq.gz"
read2="$input_dir/${sample}_R2_001_paired.fastq.gz"
if [[ -f "$read1" && -f "$read2" ]]; then
bowtie2 --threads 4 -x "$index_dir" -1 "$read1" -2 "$read2" | samtools view -bS -@ 4 - > ${outdir}/${sample}.bam
else
echo "Missing files for $sample R1 and $sample R2"
fi
done
duration=$SECONDS
echo "Process finished in $(($duration / 60)) minutes and $(($duration % 60)) seconds."
conda deactivate