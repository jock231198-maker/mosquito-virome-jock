#!/bin/bash
# ---------------------------------------------------------------------------
# star_align_worker_fast.sh — alinear al huesped con STAR y quedarse SOLO con
#                             lo no mapeado (una muestra por elemento)
#
# QUE CAMBIA RESPECTO A star_align_worker.sh
#
# 1. --outFilterMultimapNmax 10000 (el defecto de STAR es 10)
#    El genoma de Aedes aegypti es ~50% elementos transponibles: el 89.4% de las
#    lecturas cae en mas de 10 loci. Con el defecto, STAR las marca "mapped to
#    too many loci", NO las cuenta como mapeadas, y las escribe en Unmapped.out.
#    Medido en M22_S30: 20 M pares en el saco viral, frente a 1.1 M de bowtie2.
#    Subiendo el umbral pasan a "multiple loci" = huesped, que es lo correcto.
#
# 2. --outSAMmultNmax 1  <-- IMPRESCINDIBLE con lo anterior
#    STAR escribe una linea de SAM por cada locus de una lectura multimapeada.
#    Con el umbral en 10000 y el 89% multimapeando, el BAM se va a terabytes.
#    Esto limita la SALIDA a un alineamiento por lectura sin cambiar el recuento.
#
# 3. --outSAMtype None por defecto: el BAM no se quiere para nada.
#    Si resulta que suprime los Unmapped.out.mate*, el script lo detecta y avisa;
#    entonces se relanza con STAR_SAMTYPE="BAM Unsorted".
#
# 4. Comprime y renombra los no-mapeados a la convencion del pipeline.
#    STAR los escribe SIN comprimir como <sample>_Unmapped.out.mate1/2.
#
# 5. Escribe unmapped_counts/<sample>_unmapped_pairs.tsv para read_accounting.sh
#
# Variables de entorno:
#   MULTIMAP_MAX   umbral de multimapeo          (10000)
#   STAR_SAMTYPE   tipo de salida SAM/BAM        (None)
#   UNMAPPED_DIR   destino de los no-mapeados    ($SCRATCH/unmapped_fastq_star)
#
# Prep + envio:
#   source .../config.sh
#   N=$(wc -l < "$SCRATCH/samples.txt")
#   bsub -J "star[1-$N]%6" \
#        -o "$LOGS_DIR/star.%J.%I.log" -e "$LOGS_DIR/star.%J.%I.err" \
#        -q normal -n 8 -M 20000 \
#        -R "select[mem>20000] rusage[mem=20000] span[hosts=1]" \
#        "$PWD/star_align_worker_fast.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM"
#
# SALIDA: $SCRATCH/unmapped_fastq_star/<sample>_unmapped_R1.fastq.gz  (+ R2)
#         $RESULTS_DIR/star_stats/<sample>_Log.final.out
#         $RESULTS_DIR/unmapped_counts_star/<sample>_unmapped_pairs.tsv
#
# Va a carpetas propias (_star) para NO pisar los resultados de bowtie2.
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: star_align_worker_fast.sh <samples.txt> <input_dir> [genome_dir] [outbase]}"
indir="${2:?Falta input_dir}"
genome_dir="${3:-$STAR_INDEX}"
outbase="${4:-$SCRATCH/aligned}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
MULTIMAP_MAX="${MULTIMAP_MAX:-10000}"
SAM_TYPE="${STAR_SAMTYPE:-None}"          # sin comillas al usarlo: pueden ser 2 palabras
unmapdir="${UNMAPPED_DIR:-$SCRATCH/unmapped_fastq_star}"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_R1_001_paired.fastq.gz"
R2="$indir/${sample}_R2_001_paired.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: faltan R1/R2 para $sample en $indir"; exit 1; }

[[ -s "$genome_dir/SA" ]] || { echo "ERROR: indice STAR incompleto en $genome_dir"; exit 1; }

outdir="$outbase/$RESULT_FROM/${sample}"
statsdir="$RESULTS_DIR/star_stats"
countdir="$RESULTS_DIR/unmapped_counts_star"
mkdir -p "$outdir" "$statsdir" "$countdir" "$unmapdir"

activate_env "$ENV_STAR"

echo "STAR align $sample  (host: $(hostname), hilos: $THREADS)"
echo "  indice          : $genome_dir"
echo "  multimapNmax    : $MULTIMAP_MAX"
echo "  outSAMtype      : $SAM_TYPE"
echo "  no-mapeados a   : $unmapdir"

STAR --runThreadN "$THREADS" \
     --runMode alignReads \
     --genomeDir "$genome_dir" \
     --readFilesIn "$R1" "$R2" \
     --readFilesCommand zcat \
     --outFilterMultimapNmax "$MULTIMAP_MAX" \
     --outSAMmultNmax 1 \
     --outSAMtype $SAM_TYPE \
     --outReadsUnmapped Fastx \
     --outFileNamePrefix "${outdir}/${sample}_"

# --- Log final -> resultado permanente ; temporales fuera --------------------
cp "${outdir}/${sample}_Log.final.out" "$statsdir/${sample}_Log.final.out" 2>/dev/null || true
rm -rf "${outdir}/${sample}__STARtmp"

# --- Comprimir y renombrar los no-mapeados ----------------------------------
m1="${outdir}/${sample}_Unmapped.out.mate1"
m2="${outdir}/${sample}_Unmapped.out.mate2"

if [[ ! -s "$m1" || ! -s "$m2" ]]; then
  echo "ERROR: STAR no genero los Unmapped.out.mate*"
  echo "       Con --outSAMtype None puede que se supriman."
  echo "       Relanza con:  STAR_SAMTYPE='BAM Unsorted'"
  exit 1
fi

n1=$(wc -l < "$m1"); n2=$(wc -l < "$m2")
(( n1 == n2 )) || { echo "ERROR: mate1 y mate2 descuadran ($n1 vs $n2 lineas)"; exit 1; }
pairs=$(( n1 / 4 ))

if command -v pigz >/dev/null 2>&1; then COMP=(pigz -p "$THREADS" -c); else COMP=(gzip -c); fi

# .part y luego mv: un job matado no debe dejar un .gz con nombre bueno y
# contenido truncado. Es el fallo que costo relanzar el trimming.
"${COMP[@]}" "$m1" > "$unmapdir/${sample}_unmapped_R1.part"
"${COMP[@]}" "$m2" > "$unmapdir/${sample}_unmapped_R2.part"
for p in R1 R2; do
  gzip -t "$unmapdir/${sample}_unmapped_${p}.part" 2>/dev/null || {
    echo "ERROR: gzip corrupto en $p"
    rm -f "$unmapdir/${sample}_unmapped_R1.part" "$unmapdir/${sample}_unmapped_R2.part"
    exit 1; }
done
mv "$unmapdir/${sample}_unmapped_R1.part" "$unmapdir/${sample}_unmapped_R1.fastq.gz"
mv "$unmapdir/${sample}_unmapped_R2.part" "$unmapdir/${sample}_unmapped_R2.fastq.gz"
rm -f "$m1" "$m2"

printf 'sample\tunmapped_pairs\n%s\t%s\n' "$sample" "$pairs" \
  > "$countdir/${sample}_unmapped_pairs.tsv"

# --- Resumen en el log ------------------------------------------------------
L="$statsdir/${sample}_Log.final.out"
get() { sed -n "s/.*$1 |[[:space:]]*//p" "$L" | head -1; }
echo "  entrada         : $(get 'Number of input reads') pares"
echo "  unicos          : $(get 'Uniquely mapped reads %')"
echo "  multiples       : $(get '% of reads mapped to multiple loci')"
echo "  demasiados loci : $(get '% of reads mapped to too many loci')   <- debe ser ~0%"
echo "  no mapeados     : $pairs pares  -> $unmapdir/${sample}_unmapped_R1.fastq.gz"
echo "Done: $sample  ($((SECONDS / 60)) min $((SECONDS % 60)) s)"