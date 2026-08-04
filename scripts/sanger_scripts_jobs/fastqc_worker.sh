#!/bin/bash
# ---------------------------------------------------------------------------
# fastqc_worker.sh
# Procesa UN archivo fastq.gz por elemento de un job array de LSF.
#
# Cómo lanzarlo (desde un nodo de login, en tu scratch de Lustre):
#
#   cd /lustre/scratch126/tol/teams/lawniczak/users/jm79/raw_data/transcriptome/MERIDA/
#   mkdir -p logs fastqc_results
#   find raw_fastq -type f -name "*.fastq.gz" | sort > filelist.txt
#   N=$(wc -l < filelist.txt)
#
#   bsub -J "fastqc[1-$N]%10" \
#        -o "logs/fastqc.%J.%I.log" -e "logs/fastqc.%J.%I.err" \
#        -q normal -n 4 -M 2000 \
#        -R "select[mem>2000] rusage[mem=2000] span[hosts=1]" \
#        "./fastqc_worker.sh filelist.txt fastqc_results"
#
#   El %10 limita a 10 elementos corriendo a la vez (empieza pequeño).
# ---------------------------------------------------------------------------
set -euo pipefail

filelist="${1:?Uso: fastqc_worker.sh <filelist.txt> <output_dir>}"
outdir="${2:?Falta output_dir}"

# Nº de hilos: usa los CPUs que LSF te asignó (-n), por defecto 4
threads="${LSB_DJOB_NUMPROC:-4}"

# Inicializar conda y activar el entorno con FastQC instalado
module load conda
conda activate fastqc_

# Seleccionar el archivo que corresponde a este elemento del array
idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
file=$(sed -n "${idx}p" "$filelist")

echo "Array element ${idx}: ${file}  (host: $(hostname))"
mkdir -p "$outdir"

fastqc -t "$threads" "$file" -o "$outdir"

echo "FastQC terminado para: $file"
