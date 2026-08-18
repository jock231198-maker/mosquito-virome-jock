#!/bin/bash
# ---------------------------------------------------------------------------
# bowtie_map_worker_fast.sh — mapear al huesped y quedarse SOLO con lo no
#                             mapeado (una muestra por elemento)
#
# QUE CAMBIA RESPECTO A bowtie_map_worker.sh
#   El worker original manda todo el SAM a `samtools view -bS`, o sea comprime
#   ~100 GB de BAM (el 90% huesped) para luego borrarlo tras extraer los
#   no-mapeados. Su propia cabecera lo dice: "intermedio (BORRABLE)".
#
#   Aqui se usa --un-conc-gz, que hace que bowtie2 escriba directamente las
#   parejas que NO alinean de forma concordante. Se escribe ~10 GB en vez de
#   ~100, y se elimina el paso extract_unmapped entero: el nombre de salida ya
#   es el que espera `check_step.sh unmap`.
#
#   Es el mismo truco que la rama de STAR con --outReadsUnmapped Fastx.
#
# MATIZ IMPORTANTE
#   --un-conc-gz guarda las parejas que no alinean CONCORDANTEMENTE, que son
#   algo mas que aquellas con ambos extremos sin mapear (el criterio del
#   extract_unmapped clasico, samtools -f 12). En las catas la diferencia era
#   ~1.4%. Para virome peca por el lado seguro: es preferible arrastrar algo de
#   huesped que perder una read viral cuya pareja mapeo por casualidad.
#
# NO NECESITA samtools.
#
# Prep + envio:
#   source .../config.sh
#   N=$(wc -l < "$SCRATCH/samples.txt")
#   bsub -J "btmap[1-$N]%10" \
#        -o "$LOGS_DIR/btmap.%J.%I.log" -e "$LOGS_DIR/btmap.%J.%I.err" \
#        -q normal -n 16 -M 12000 \
#        -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" \
#        "$PWD/bowtie_map_worker_fast.sh $SCRATCH/samples.txt $SCRATCH/trimmed/$RESULT_FROM"
#
# SALIDA (intermedia): $SCRATCH/unmapped_fastq/<sample>_unmapped_R1.fastq.gz
#                      $SCRATCH/unmapped_fastq/<sample>_unmapped_R2.fastq.gz
# SE CONSERVA        : $RESULTS_DIR/mapping_stats/<sample>_bowtie2_summary.txt
#                      $RESULTS_DIR/unmapped_counts/<sample>_unmapped_pairs.tsv
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: bowtie_map_worker_fast.sh <samples.txt> <input_dir> [index] [outdir]}"
indir="${2:?Falta input_dir}"
index="${3:-$BT2_INDEX}"
outdir="${4:-$SCRATCH/unmapped_fastq}"

THREADS="${LSB_DJOB_NUMPROC:-8}"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_R1_001_paired.fastq.gz"
R2="$indir/${sample}_R2_001_paired.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: faltan R1/R2 para $sample en $indir"; exit 1; }

statsdir="$RESULTS_DIR/mapping_stats"
countdir="$RESULTS_DIR/unmapped_counts"

# Staging aparte: un job matado no debe dejar un .gz con nombre bueno y
# contenido truncado. Es el fallo que nos costo relanzar el trimming.
# Va FUERA de $outdir porque check_step.sh busca con find recursivo.
stage="${BTMAP_STAGE:-$SCRATCH/.btmap_stage}"
mkdir -p "$outdir" "$statsdir" "$countdir" "$stage"

O1="$outdir/${sample}_unmapped_R1.fastq.gz"
O2="$outdir/${sample}_unmapped_R2.fastq.gz"
T1="$stage/${sample}_unmapped_R1.fastq.gz"
T2="$stage/${sample}_unmapped_R2.fastq.gz"
rm -f "$T1" "$T2"

activate_env "$ENV_BOWTIE2"

echo "Mapping $sample  (host: $(hostname), hilos: $THREADS)"
echo "  indice : $index"
echo "  salida : $outdir"

# El resumen de bowtie2 va a stderr -> se guarda como resultado permanente.
# El SAM se descarta: solo interesan las parejas no concordantes.
bowtie2 --threads "$THREADS" \
        -x "$index" \
        -1 "$R1" -2 "$R2" \
        --un-conc-gz "$stage/${sample}_unmapped_R%.fastq.gz" \
        -S /dev/null \
        2> "$statsdir/${sample}_bowtie2_summary.txt"

# --- Verificacion antes de dar por buenos los ficheros ----------------------
for f in "$T1" "$T2"; do
  [[ -s "$f" ]] || { echo "ERROR: no se genero $f"; rm -f "$T1" "$T2"; exit 1; }
  if ! gzip -t "$f" 2>/dev/null; then
    echo "ERROR: $f no pasa gzip -t"; rm -f "$T1" "$T2"; exit 1
  fi
done

# Mismo numero de reads en R1 y R2, o el emparejamiento esta roto
n1=$(gzip -cd "$T1" | wc -l); n2=$(gzip -cd "$T2" | wc -l)
if (( n1 != n2 )); then
  echo "ERROR: R1 tiene $((n1/4)) reads y R2 $((n2/4)) -> desemparejado"
  rm -f "$T1" "$T2"; exit 1
fi
pairs=$((n1 / 4))

mv "$T1" "$O1"
mv "$T2" "$O2"

# --- Conteo para read_accounting.sh -----------------------------------------
# La columna unmapped_pairs sale de aqui. Nadie mas escribia este fichero.
printf 'sample\tunmapped_pairs\n%s\t%s\n' "$sample" "$pairs" \
  > "$countdir/${sample}_unmapped_pairs.tsv"

# --- Resumen en el log ------------------------------------------------------
rate=$(sed -n 's/^\([0-9.]*\)% overall alignment rate.*/\1/p' "$statsdir/${sample}_bowtie2_summary.txt" | head -1)
total=$(sed -n 's/^\([0-9]*\) reads; of these:.*/\1/p'          "$statsdir/${sample}_bowtie2_summary.txt" | head -1)
echo "  pares procesados : ${total:-?}"
echo "  mapean al huesped: ${rate:-?}%"
echo "  no concordantes  : $pairs pares  -> $O1"
echo "Done: $sample  ($((SECONDS / 60)) min $((SECONDS % 60)) s)"
