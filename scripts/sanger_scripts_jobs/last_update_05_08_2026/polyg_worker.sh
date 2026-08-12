#!/bin/bash
# ---------------------------------------------------------------------------
# polyg_worker.sh — elimina las colas polyG (una muestra por elemento)
#
# POR QUE EXISTE ESTE PASO
#   La quimica de dos colores (NovaSeq/NextSeq) llama G cuando no hay senal.
#   En este lote el polyG es la secuencia sobrerrepresentada nº1: aparece en
#   las 23 muestras, 23.2 M ocurrencias, 1.23% de todas las reads, con colas
#   de hasta el 12.8% en M74_S21_R2.
#
#   Trimmomatic NO puede quitarlo: no tiene recorte de polyG, y el recorte por
#   calidad tampoco sirve porque esas G suelen venir con Q alto. Si no se trata,
#   las colas no mapean al huesped, caen en el saco de "no mapeadas", entran a
#   SPAdes y salen como contigs de baja complejidad.
#
#   fastp aqui hace UNA sola cosa. Las tres --disable_* apagan todo lo demas
#   para no pisar el trabajo de Trimmomatic ni cambiar sus parametros.
#
# Prep + envio:
#   source .../config.sh
#   ls "$SCRATCH/cat_fastq"/*_R1_001.fastq.gz | xargs -n1 basename | cut -d_ -f1-2 \
#       | sort -u | grep -vx 'M41_S4' > "$SCRATCH/samples.txt"
#   N=$(wc -l < "$SCRATCH/samples.txt")
#
#   bsub -J "polyg[1-$N]%10" \
#        -o "$LOGS_DIR/polyg.%J.%I.log" -e "$LOGS_DIR/polyg.%J.%I.err" \
#        -q normal -n 4 -M 4000 \
#        -R "select[mem>4000] rusage[mem=4000] span[hosts=1]" \
#        "$PWD/polyg_worker.sh $SCRATCH/samples.txt $SCRATCH/cat_fastq"
#
# SALIDA (intermedia, en SCRATCH): $SCRATCH/nopolyg/
#   <sample>_R1_001.fastq.gz  /  <sample>_R2_001.fastq.gz   <- entrada del trimming
# INFORMES (se conservan): $RESULTS_DIR/polyg_reports/<sample>.json + .html
#   El .json lo agrega MultiQC: da el antes/despues de cada muestra.
#
# Cadena completa:  cat_fastq -> nopolyg -> trimmed
# Se puede borrar nopolyg/ en cuanto el trimming termine (es regenerable).
# ---------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: polyg_worker.sh <samples.txt> <input_dir> [output_dir]}"
indir="${2:?Falta input_dir}"
outdir="${3:-$SCRATCH/nopolyg}"

# Longitud minima de la carrera de G para considerarla cola. fastp usa 10 por
# defecto; 8 es algo mas agresivo y aqui compensa, porque las colas cortas son
# las que mas abundan.
POLYG_MIN="${POLYG_MIN:-8}"

# fastp no escala mas alla de ~16 hilos
THREADS="${LSB_DJOB_NUMPROC:-4}"
(( THREADS > 16 )) && THREADS=16

activate_env "$ENV_FASTP"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_R1_001.fastq.gz"
R2="$indir/${sample}_R2_001.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: faltan R1/R2 para $sample en $indir"; exit 1; }

repdir="$RESULTS_DIR/polyg_reports"
mkdir -p "$outdir" "$repdir"

O1="$outdir/${sample}_R1_001.fastq.gz"
O2="$outdir/${sample}_R2_001.fastq.gz"

echo "polyG $sample  (host: $(hostname), hilos: $THREADS, poly_g_min_len: $POLYG_MIN)"

# Se escribe a .tmp y se renombra al final: un job cortado a medias nunca deja
# un .gz truncado con pinta de bueno. Mismo criterio que compress.sh.
fastp \
    -i "$R1" -I "$R2" \
    -o "$O1.tmp" -O "$O2.tmp" \
    --trim_poly_g --poly_g_min_len "$POLYG_MIN" \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    --thread "$THREADS" \
    --compression 6 \
    --json "$repdir/${sample}.json" \
    --html "$repdir/${sample}.html" \
    --report_title "polyG $sample"

# Verificacion: que los dos .gz esten integros antes de darlos por buenos
for f in "$O1.tmp" "$O2.tmp"; do
  if ! gzip -t "$f" 2>/dev/null; then
    echo "ERROR: $f no pasa gzip -t"; rm -f "$O1.tmp" "$O2.tmp"; exit 1
  fi
done

mv "$O1.tmp" "$O1"
mv "$O2.tmp" "$O2"

# Resumen legible en el log, por si no quieres abrir el JSON
if command -v python3 >/dev/null 2>&1; then
  python3 - "$repdir/${sample}.json" <<'PY' || true
import json, sys
d = json.load(open(sys.argv[1]))
b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
pg = d.get("polyg_trimming", {}).get("total_polyg_trimmed_reads", "?")
print(f"  reads   : {b['total_reads']:,} -> {a['total_reads']:,}")
print(f"  bases   : {b['total_bases']:,} -> {a['total_bases']:,} "
      f"({100*(b['total_bases']-a['total_bases'])/b['total_bases']:.2f}% recortado)")
print(f"  long.med: {b['read1_mean_length']} -> {a['read1_mean_length']} (R1)")
print(f"  reads con polyG recortado: {pg}")
PY
fi

echo "Done: $sample"
