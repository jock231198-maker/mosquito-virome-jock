#!/bin/bash
# ---------------------------------------------------------------------------
# assembly_trinity.sh — ensamblaje con Trinity (una muestra por elemento)
#
# Mismo contrato de salida que assembly_spades.sh:
#   $SCRATCH/trinity/<sample>/contigs.fasta          <- quast, read_accounting
#   $SCRATCH/trinity/<sample>/contigs_1000bp.fasta   <- diamond
#   $SCRATCH/trinity/<sample>/<sample>_final.fasta   <- genomad, checkv
#   $SCRATCH/trinity/<sample>/Trinity.timing, gene_trans_map
#   $RESULTS_DIR/assembly_stats_trinity/<sample>_assembly.tsv
#
# CUATRO COSAS QUE HAY QUE SABER ANTES DE LANZAR ESTO
#
# 1. INODES. La fase Chrysalis parte el grafo en decenas de miles de ficheros
#    diminutos. En Lustre con cuota de inodes eso tumba el trabajo -- y de paso
#    a los demas. Por eso TODO ocurre en el disco local del nodo y a Lustre solo
#    llega el resultado. No quites esto por comodidad.
#
# 2. --full_cleanup. Trinity borra su directorio de trabajo y deja solo el FASTA
#    final. Es la otra mitad de la defensa contra el punto 1. El precio es que
#    no se puede reanudar un trabajo interrumpido: hay que repetirlo entero.
#    Ponlo a 0 solo si estas depurando una muestra concreta.
#
# 3. EL NOMBRE DEL DIRECTORIO DE SALIDA DEBE CONTENER "trinity". Trinity aborta
#    si no. Es una comprobacion suya, no un capricho de este script.
#
# 4. NO HAY COBERTURA POR CONTIG. Las cabeceras de Trinity son
#       >TRINITY_DN1000_c0_g1_i1 len=329 path=[...]
#    -- longitud si, cobertura no. SPAdes y MEGAHIT si la dan. Si la necesitas
#    para separar virus abundante de fragmento espurio, hay que mapear las
#    lecturas de vuelta (align_and_estimate_abundance.pl, de la propia Trinity).
#    La columna 'coverage' de contigs_table.tsv sale vacia a proposito.
#
# Variables de entorno:
#   TRINITY_MEM       --max_memory en GB                       (32)
#   TRINITY_MINLEN    --min_contig_length                      (200)
#   TRINITY_NORM      1 normalizacion in silico (defecto de Trinity), 0 la apaga (1)
#   TRINITY_SS        "" | RF | FR   protocolo de hebra        ("")
#   FULL_CLEANUP      1 borra el directorio de trabajo          (1)
#   MIN_LEN           longitud del filtro posterior            (1000)
#
# OJO: TRINITY_MEM por debajo del -M del bsub.
#
# Prep + envio (SMOKE de dos muestras, una pequena y una grande):
#   source .../config.sh
#   printf '%s\n' M36_S17 M11_S3 > "$SCRATCH/trin_smoke.txt"
#   bsub -J "trinsmoke[1-2]" \
#        -o "$LOGS_DIR/trinsmoke.%J.%I.log" -e "$LOGS_DIR/trinsmoke.%J.%I.err" \
#        -q long -n 8 -M 40000 \
#        -R "select[mem>40000] rusage[mem=40000] span[hosts=1]" \
#        "$SCRIPTS_DIR/assembly_trinity.sh $SCRATCH/trin_smoke.txt $SCRATCH/unmapped_fastq"
#
# La cola es 'long' y no 'normal' a proposito: Trinity son horas, no minutos.
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

samples="${1:?Uso: assembly_trinity.sh <samples.txt> <input_dir> [outbase]}"
indir="${2:?Falta input_dir}"
outbase="${3:-$SCRATCH/trinity}"

THREADS="${LSB_DJOB_NUMPROC:-8}"
TRINITY_MEM="${TRINITY_MEM:-32}"          # GB
TRINITY_MINLEN="${TRINITY_MINLEN:-200}"
TRINITY_NORM="${TRINITY_NORM:-1}"
TRINITY_SS="${TRINITY_SS:-}"
FULL_CLEANUP="${FULL_CLEANUP:-1}"
MIN_LEN="${MIN_LEN:-1000}"

idx="${LSB_JOBINDEX:?Este script debe enviarse como job array}"
sample=$(sed -n "${idx}p" "$samples")
[[ -n "$sample" ]] || { echo "ERROR: linea $idx vacia en $samples"; exit 1; }

R1="$indir/${sample}_unmapped_R1.fastq.gz"
R2="$indir/${sample}_unmapped_R2.fastq.gz"
[[ -f "$R1" && -f "$R2" ]] || { echo "ERROR: faltan R1/R2 para $sample en $indir"; exit 1; }

final_dir="$outbase/$sample"
statsdir="${ASM_STATS_DIR:-$RESULTS_DIR/assembly_stats_trinity}"
mkdir -p "$outbase" "$statsdir"

# --- Idempotencia -----------------------------------------------------------
if [[ -s "$final_dir/contigs.fasta" ]]; then
  echo "Ya existe $final_dir/contigs.fasta ($(grep -c '^>' "$final_dir/contigs.fasta" || true) contigs). No se rehace."
  echo "Para rehacerlo:  rm -rf $final_dir"
  exit 0
fi

# --- Trabajo en disco local del nodo ----------------------------------------
work="${TMPDIR:-/tmp}/trin_${sample}_${LSB_JOBID:-$$}"
rm -rf "$work"; mkdir -p "$work"
cleanup() { rm -rf "$work"; }
trap cleanup EXIT

# Espacio libre en el disco local. Trinity necesita holgura y es mejor enterarse
# ahora que a las seis horas de ejecucion.
avail_gb=$(df -BG --output=avail "$work" 2>/dev/null | tail -1 | tr -dc '0-9')
if [[ -n "${avail_gb:-}" ]] && (( avail_gb < 50 )); then
  echo "AVISO: solo ${avail_gb} GB libres en $(dirname "$work"). Trinity puede quedarse sin sitio."
fi

# Requisito de Trinity: el nombre debe contener "trinity".
trin_out="$work/trinity_${sample}"

activate_env "$ENV_TRINITY"
command -v Trinity >/dev/null || { echo "ERROR: Trinity no esta en el PATH (ojo: T mayuscula)"; exit 1; }

# --- Opciones ---------------------------------------------------------------
opts=()
[[ "$TRINITY_NORM"  == "0" ]] && opts+=( --no_normalize_reads )
[[ "$FULL_CLEANUP"  == "1" ]] && opts+=( --full_cleanup )
[[ -n "$TRINITY_SS" ]]        && opts+=( --SS_lib_type "$TRINITY_SS" )

echo "==========================================================="
echo "Trinity $sample   (host: $(hostname))"
echo "  hilos     : $THREADS      memoria: ${TRINITY_MEM}G"
echo "  normaliz. : $([[ "$TRINITY_NORM" == "1" ]] && echo "si (defecto, max_cov 200)" || echo "NO")"
echo "  hebra     : ${TRINITY_SS:-<sin hebra>}"
echo "  cleanup   : $FULL_CLEANUP"
echo "  entrada   : $R1"
echo "  trabajo   : $work   (disco local, ${avail_gb:-?} GB libres)"
Trinity --version 2>&1 | head -2 || true
echo "==========================================================="

Trinity \
    --seqType fq \
    --left  "$R1" \
    --right "$R2" \
    --CPU "$THREADS" \
    --max_memory "${TRINITY_MEM}G" \
    --min_contig_length "$TRINITY_MINLEN" \
    --output "$trin_out" \
    ${opts[@]+"${opts[@]}"}

# --- Localizar el FASTA -----------------------------------------------------
# Trinity escribe el resultado JUNTO al directorio, no dentro:
#   con --full_cleanup   -> trinity_<sample>.Trinity.fasta  (y borra el directorio)
#   sin --full_cleanup   -> trinity_<sample>/Trinity.fasta
# El nombre ha cambiado entre versiones, asi que se prueban las dos formas.
asm=""
for cand in "${trin_out}.Trinity.fasta" "$trin_out/Trinity.fasta"; do
  [[ -s "$cand" ]] && { asm="$cand"; break; }
done
[[ -n "$asm" ]] || {
  echo "ERROR: no aparece el FASTA de Trinity. Contenido de $work:"
  ls -la "$work" || true; exit 1; }
echo "  ensamblaje encontrado: $asm"

stage="$work/publish"
mkdir -p "$stage"
cp "$asm" "$stage/contigs.fasta"
for extra in "${asm}.gene_trans_map" "${trin_out}/Trinity.timing" "${trin_out}.Trinity.fasta.gene_trans_map"; do
  [[ -f "$extra" ]] && cp "$extra" "$stage/" || true
done

n_all=$(grep -c '^>' "$stage/contigs.fasta" || true)
(( n_all > 0 )) || { echo "ERROR: contigs.fasta vacio"; exit 1; }

# --- Filtro por longitud ----------------------------------------------------
filt="$stage/contigs_${MIN_LEN}bp.fasta"
if command -v seqkit >/dev/null 2>&1; then
  seqkit seq -m "$MIN_LEN" "$stage/contigs.fasta" > "$filt"
else
  awk -v L="$MIN_LEN" '
    /^>/ { if (h != "" && length(s) >= L) print h ORS s; h=$0; s=""; next }
         { s = s $0 }
    END  { if (h != "" && length(s) >= L) print h ORS s }' \
    "$stage/contigs.fasta" > "$filt"
fi
n_filt=$(grep -c '^>' "$filt" || true)

cp "$filt" "$stage/${sample}_final.fasta"

# --- Tabla de contigs, formato Trinity --------------------------------------
# >TRINITY_DN1000_c0_g1_i1 len=329 path=[...]
# La columna coverage va vacia: Trinity no la reporta (ver punto 4 de la cabecera).
# 'isoform' permite agrupar: _i1, _i2... del mismo _g son variantes del mismo gen,
# algo que ni SPAdes ni MEGAHIT producen y que infla el recuento de contigs.
{
  printf 'sample\tnode\tlength_bp\tcoverage\tgene\n'
  awk -v s="$sample" '/^>/ {
      id = substr($1, 2)
      len = ""
      for (i = 2; i <= NF; i++) if ($i ~ /^len=/) len = substr($i, 5)
      gene = id
      sub(/_i[0-9]+$/, "", gene)
      print s "\t" id "\t" len "\t\t" gene
    }' "$stage/contigs.fasta"
} > "$stage/contigs_table.tsv"

n_genes=$(tail -n +2 "$stage/contigs_table.tsv" | cut -f5 | sort -u | wc -l)

# --- Publicar de forma atomica ----------------------------------------------
rm -rf "$final_dir"
mkdir -p "$(dirname "$final_dir")"
mv "$stage" "$final_dir"

# --- Estadisticas -----------------------------------------------------------
longest=$(awk '/^>/{if(l>m)m=l; l=0; next}{l+=length($0)}END{if(l>m)m=l; print m+0}' \
          "$final_dir/contigs.fasta")
total_bp=$(awk '!/^>/{t+=length($0)}END{print t+0}' "$final_dir/contigs.fasta")

printf 'sample\tassembler\tmode\tcontigs_total\tcontigs_ge%s\ttotal_bp\tlongest_bp\tgenes\n' "$MIN_LEN" \
  > "$statsdir/${sample}_assembly.tsv"
printf '%s\ttrinity\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "$sample" "$([[ "$TRINITY_NORM" == "1" ]] && echo norm || echo nonorm)" \
  "$n_all" "$n_filt" "$total_bp" "$longest" "$n_genes" \
  >> "$statsdir/${sample}_assembly.tsv"

echo
echo "  transcritos totales : $n_all   (en $n_genes 'genes')"
echo "  >=${MIN_LEN} bp          : $n_filt"
echo "  bases totales       : $total_bp"
echo "  mas largo           : $longest bp"
echo "  salida              : $final_dir"
echo "Done: $sample  ($((SECONDS/3600))h $((SECONDS%3600/60))m $((SECONDS%60))s)"
