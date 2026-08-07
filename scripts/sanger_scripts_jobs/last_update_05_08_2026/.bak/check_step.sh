#!/bin/bash
# ---------------------------------------------------------------------------
# check_step.sh — verifica que un paso del pipeline termino COMPLETO.
#
# El fallo mas caro en un job array no es el que peta con error, es el que
# desaparece en silencio: 3 elementos de 200 mueren por memoria, nadie lo ve,
# y sigues el pipeline con 197 muestras creyendo que tienes 200.
#
# Comprueba tres cosas:
#   1. Los logs de LSF: hay Exited / TERM_MEMLIMIT / TERM_RUNLIMIT?
#   2. Cuenta de entradas vs cuenta de salidas
#   3. Salidas vacias o sospechosamente pequenas
#
# Uso:
#   ./check_step.sh fastqc
#   ./check_step.sh trim
#   ./check_step.sh btmap
#   ./check_step.sh --list
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

step="${1:-}"

if [[ -z "$step" || "$step" == "--list" ]]; then
  cat <<'EOF'
Uso: ./check_step.sh <paso>

Pasos:
  fastqc    FastQC sobre crudos
  trim      Trimmomatic
  btmap     Mapeo con bowtie2
  unmap     Extraccion de no-mapeados
  star      Alineamiento STAR
  spades    Ensamblaje
  quast     QC de ensamblaje
  genomad   geNomad
  checkv    CheckV
  diamond   DIAMOND blastx
EOF
  exit 0
fi

ERRORS=0; WARNS=0
err()  { echo "  [ERROR] $*"; ERRORS=$((ERRORS+1)); }
warn() { echo "  [warn ] $*"; WARNS=$((WARNS+1)); }
ok()   { echo "  [ok   ] $*"; }

# --- Configuracion por paso -------------------------------------------------
# in_key/out_key indican de donde sale el nombre de muestra:
#   file = del nombre del fichero      (M1_S1.bam        -> M1_S1)
#   dir  = del nombre de su carpeta    (spades/M1_S1/contigs.fasta -> M1_S1)
case "$step" in
  fastqc)
    logpre="fastqc"
    in_dir="$DATA_DIR";                      in_pat="*_R1_*.fastq.gz";        in_key=file
    out_dir="$RESULTS_DIR/fastqc_raw";       out_pat="*_R1*_fastqc.zip";      out_key=file
    min_bytes=100000 ;;
  trim)
    logpre="trim"
    in_dir="$DATA_DIR";                      in_pat="*_R1_*.fastq.gz";        in_key=file
    out_dir="$SCRATCH/trimmed/$RESULT_FROM"; out_pat="*_R1_001_paired.fastq.gz"; out_key=file
    min_bytes=1000000 ;;
  btmap)
    logpre="btmap"
    in_dir="$SCRATCH/trimmed/$RESULT_FROM";  in_pat="*_R1_001_paired.fastq.gz"; in_key=file
    out_dir="$SCRATCH/mapped/$RESULT_FROM";  out_pat="*.bam";                 out_key=file
    min_bytes=1000000 ;;
  unmap)
    logpre="unmap"
    in_dir="$SCRATCH/mapped";                in_pat="*.bam";                  in_key=file
    out_dir="$SCRATCH/unmapped_fastq";       out_pat="*_unmapped_R1.fastq.gz"; out_key=file
    min_bytes=1000 ;;
  star)
    logpre="star"
    in_dir="$SCRATCH/trimmed/$RESULT_FROM";  in_pat="*_R1_001_paired.fastq.gz"; in_key=file
    out_dir="$SCRATCH/aligned";              out_pat="*_Aligned.out.bam";     out_key=file
    min_bytes=1000000 ;;
  spades)
    logpre="spades"
    in_dir="$SCRATCH/unmapped_fastq";        in_pat="*_unmapped_R1.fastq.gz"; in_key=file
    out_dir="$SCRATCH/spades";               out_pat="contigs.fasta";         out_key=dir
    min_bytes=1000 ;;
  quast)
    logpre="quast"
    in_dir="$SCRATCH/spades";                in_pat="contigs.fasta";          in_key=dir
    out_dir="$RESULTS_DIR/quast";            out_pat="report.tsv";            out_key=dir
    min_bytes=100 ;;
  genomad)
    logpre="genomad"
    in_dir="$SCRATCH/spades";                in_pat="*_final.fasta";          in_key=file
    out_dir="$RESULTS_DIR/genomad";          out_pat="*_summary";             out_key=dir
    min_bytes=1 ;;
  checkv)
    logpre="checkv"
    in_dir="$SCRATCH/spades";                in_pat="*_final.fasta";          in_key=file
    out_dir="$RESULTS_DIR/checkv";           out_pat="quality_summary.tsv";   out_key=dir
    min_bytes=100 ;;
  diamond)
    logpre="dmnd"
    in_dir="$SCRATCH/spades";                in_pat="contigs_1000bp.fasta";   in_key=dir
    out_dir="$RESULTS_DIR/diamond";          out_pat="*_diamond_extended.tsv"; out_key=file
    min_bytes=1 ;;
  *)
    echo "Paso desconocido: $step  (usa --list)"; exit 1 ;;
esac

# Extrae el nombre de muestra de una ruta.
#   file: basename -> quita desde el primer punto -> primeros 2 campos por '_'
#   dir : nombre de la carpeta que lo contiene
sample_key() {
  local path="$1" mode="$2" b
  if [[ "$mode" == "dir" ]]; then
    basename "$(dirname "$path")"
  else
    b=$(basename "$path")
    b="${b%%.*}"                 # M1_S1.bam -> M1_S1 ; M1_S1_R1_001_paired.fastq.gz -> M1_S1_R1_001_paired
    echo "$b" | cut -d_ -f1-2    # -> M1_S1
  fi
}

keys_of() {   # keys_of <dir> <pat> <mode>
  find "$1" -name "$2" 2>/dev/null | while IFS= read -r p; do
    sample_key "$p" "$3"
  done | sort -u
}

echo "==========================================================="
echo "check_step.sh  paso: $step     $(date)"
echo "==========================================================="
echo

# --- 1. Logs de LSF ---------------------------------------------------------
echo "== 1. Logs de LSF ($LOGS_DIR/$logpre.*)"
shopt -s nullglob
logs=( "$LOGS_DIR/$logpre".*.log )
shopt -u nullglob

if [[ ${#logs[@]} -eq 0 ]]; then
  warn "no hay logs con prefijo '$logpre' (paso no lanzado todavia?)"
else
  echo "  ${#logs[@]} logs encontrados"
  # LSF escribe "Successfully completed." o "Exited with exit code N"
  n_ok=$(grep -l "Successfully completed" "${logs[@]}" 2>/dev/null | wc -l)
  n_exit=$(grep -l "Exited with" "${logs[@]}" 2>/dev/null | wc -l)
  n_mem=$(grep -l "TERM_MEMLIMIT" "${logs[@]}" 2>/dev/null | wc -l)
  n_run=$(grep -l "TERM_RUNLIMIT" "${logs[@]}" 2>/dev/null | wc -l)

  echo "  completados OK : $n_ok"
  [[ $n_exit -gt 0 ]] && err "salieron con error: $n_exit"
  [[ $n_mem  -gt 0 ]] && err "matados por MEMORIA (TERM_MEMLIMIT): $n_mem  -> sube el -M"
  [[ $n_run  -gt 0 ]] && err "matados por TIEMPO (TERM_RUNLIMIT): $n_run  -> cola mas larga"

  if [[ $n_exit -gt 0 || $n_mem -gt 0 ]]; then
    echo "  elementos fallidos:"
    grep -l -e "Exited with" -e "TERM_MEMLIMIT" "${logs[@]}" 2>/dev/null \
      | sed 's/^/         /' | head -10
  fi
fi
echo

# --- 2. Entradas vs salidas -------------------------------------------------
echo "== 2. Muestras con entrada vs muestras con salida"
mapfile -t in_keys  < <(keys_of "$in_dir"  "$in_pat"  "$in_key")
mapfile -t out_keys < <(keys_of "$out_dir" "$out_pat" "$out_key")
n_in=${#in_keys[@]}; n_out=${#out_keys[@]}

echo "  entradas ($in_dir / $in_pat):  $n_in muestras"
echo "  salidas  ($out_dir / $out_pat): $n_out muestras"

if [[ $n_in -eq 0 ]]; then
  warn "no hay entradas: revisa la ruta"
else
  missing=$(comm -23 \
      <(printf '%s\n' "${in_keys[@]}") \
      <(printf '%s\n' ${out_keys[@]+"${out_keys[@]}"}) )
  extra=$(comm -13 \
      <(printf '%s\n' "${in_keys[@]}") \
      <(printf '%s\n' ${out_keys[@]+"${out_keys[@]}"}) )

  if [[ -n "$missing" ]]; then
    err "$(printf '%s\n' "$missing" | wc -l) muestras SIN salida:"
    printf '%s\n' "$missing" | sed 's/^/         /' | head -20
  else
    ok "todas las $n_in muestras tienen salida"
  fi
  if [[ -n "$extra" ]]; then
    warn "salidas sin entrada correspondiente (restos de otra corrida?):"
    printf '%s\n' "$extra" | sed 's/^/         /' | head -10
  fi
fi
echo

# --- 3. Salidas vacias o minusculas -----------------------------------------
echo "== 3. Salidas vacias o sospechosas (< $(numfmt --to=iec-i --suffix=B "$min_bytes" 2>/dev/null || echo "${min_bytes}B"))"
n_small=0
while IFS= read -r f; do
  [[ -z "$f" ]] && continue
  sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
  if [[ $sz -lt $min_bytes ]]; then
    warn "$(basename "$(dirname "$f")")/$(basename "$f"): $(numfmt --to=iec-i --suffix=B "$sz" 2>/dev/null || echo "${sz}B")"
    n_small=$((n_small+1))
  fi
done < <(find "$out_dir" -name "$out_pat" -type f 2>/dev/null)
[[ $n_small -eq 0 ]] && ok "ninguna salida vacia o minuscula"
echo

# --- 4. Consumo real de memoria (ayuda a afinar el -M) ----------------------
echo "== 4. Memoria usada realmente"
if [[ ${#logs[@]} -gt 0 ]]; then
  # LSF escribe "Max Memory : N MB" en el resumen del log
  grep -h "Max Memory" "${logs[@]}" 2>/dev/null \
    | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) {print $i; break}}' \
    | sort -n | awk '
        {v[NR]=$1}
        END{
          if(NR==0){print "  (LSF no reporto Max Memory)"; exit}
          printf "  min %s MB | mediana %s MB | max %s MB  (n=%d)\n", v[1], v[int((NR+1)/2)], v[NR], NR
          printf "  sugerencia -M: %d  (max x1.3)\n", v[NR]*1.3
        }'
else
  echo "  (sin logs)"
fi
echo

echo "==========================================================="
echo "RESUMEN paso '$step': $ERRORS errores, $WARNS avisos"
if [[ $ERRORS -gt 0 ]]; then
  echo "Relanza SOLO los elementos que fallaron antes de continuar."
  echo "Truco: bsub -J \"${logpre}rerun[i,j,k]\" ... (indices concretos entre corchetes)"
elif [[ ${#logs[@]} -eq 0 && $n_out -eq 0 ]]; then
  echo "Este paso NO se ha ejecutado todavia (sin logs y sin salidas)."
elif [[ $WARNS -gt 0 ]]; then
  echo "Sin errores, pero revisa los avisos de arriba antes de continuar."
else
  echo "Paso completo. Puedes continuar."
fi
echo "==========================================================="

exit $(( ERRORS > 0 ? 1 : 0 ))
