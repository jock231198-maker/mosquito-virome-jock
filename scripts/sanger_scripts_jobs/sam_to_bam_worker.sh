#!/bin/bash
# ---------------------------------------------------------------------------
# sam_to_bam_worker.sh  — convertir SAM -> BAM (un archivo por elemento)
#
# Solo lo necesitas si algún paso deja .sam sueltos. (bowtie_map_worker ya
# entrega BAM, y star_align_worker usa --outSAMtype BAM, así que normalmente
# no hace falta.)
#
# Prep + envío:
#   find aligned -type f -name "*.sam" | sort > sam_list.txt
#   N=$(wc -l < sam_list.txt)
#   bsub -J "sam2bam[1-$N]%10" -o "logs/sam2bam.%J.%I.log" -e "logs/sam2bam.%J.%I.err" \
#        -q normal -n 2 -M 4000 \
#        -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "./sam_to_bam_worker.sh sam_list.txt bam_out"
# ---------------------------------------------------------------------------
set -euo pipefail

sam_list="${1:?Uso: sam_to_bam_worker.sh <sam_list.txt> <output_dir>}"
outdir="${2:?Falta output_dir}"

THREADS="${LSB_DJOB_NUMPROC:-2}"
CONDA_ENV="${CONDA_ENV:-bowtie2}"   # cualquier env que tenga samtools

eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

idx="${LSB_JOBINDEX:?Enviar como job array}"
sam=$(sed -n "${idx}p" "$sam_list")
sample=$(basename "$sam" .sam)
mkdir -p "$outdir"

echo "SAM->BAM $sample (host: $(hostname))"
# Lee del .sam ENCONTRADO (tu versión leía por error desde la carpeta de salida)
samtools view -bS -@ "$THREADS" "$sam" > "$outdir/${sample}.bam"
echo "Done: $outdir/${sample}.bam"
