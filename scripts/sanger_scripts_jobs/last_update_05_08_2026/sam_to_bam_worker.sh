#!/bin/bash
# ---------------------------------------------------------------------------
# sam_to_bam_worker.sh — convertir SAM -> BAM (un archivo por elemento)
#
# Solo lo necesitas si algún paso deja .sam sueltos. (bowtie_map_worker ya
# entrega BAM y star_align_worker usa --outSAMtype BAM, así que normalmente
# no hace falta.)
#
# Prep + envío:
#   source .../config.sh
#   find "$SCRATCH/aligned" -type f -name "*.sam" | sort > "$SCRATCH/sam_list.txt"
#   N=$(wc -l < "$SCRATCH/sam_list.txt")
#   bsub -J "sam2bam[1-$N]%10" \
#        -o "$LOGS_DIR/sam2bam.%J.%I.log" -e "$LOGS_DIR/sam2bam.%J.%I.err" \
#        -q normal -n 2 -M 4000 \
#        -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "./sam_to_bam_worker.sh $SCRATCH/sam_list.txt"
#
# El .sam original queda BORRABLE una vez validado el .bam (ver cleanup.sh).
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

sam_list="${1:?Uso: sam_to_bam_worker.sh <sam_list.txt> [output_dir]}"
outdir="${2:-$SCRATCH/bam_out}"

THREADS="${LSB_DJOB_NUMPROC:-2}"
activate_env "$ENV_SAMTOOLS"

idx="${LSB_JOBINDEX:?Enviar como job array}"
sam=$(sed -n "${idx}p" "$sam_list")
[[ -n "$sam" ]] || { echo "ERROR: linea $idx vacia en $sam_list"; exit 1; }
sample=$(basename "$sam" .sam)
mkdir -p "$outdir"

echo "SAM->BAM $sample (host: $(hostname))"
# Lee del .sam ENCONTRADO (tu versión original leía por error desde la carpeta de salida)
samtools view -bS -@ "$THREADS" "$sam" > "$outdir/${sample}.bam"
echo "Done: $outdir/${sample}.bam"
