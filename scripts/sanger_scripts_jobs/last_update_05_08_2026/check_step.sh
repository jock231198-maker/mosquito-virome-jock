#!/bin/bash
# ---------------------------------------------------------------------------
# check_step.sh — verifica que un paso del pipeline termino COMPLETO.
#
# El fallo mas caro en un job array no es el que peta con error, es el que
# desaparece en silencio: 3 elementos de 46 mueren por memoria, nadie lo ve,
# y sigues el pipeline con 20 muestras creyendo que tienes 23.
#
# Comprueba cinco cosas:
#   1. Logs de LSF: Exited / TERM_MEMLIMIT / TERM_RUNLIMIT
#   2. Cuantos elementos se lanzaron vs cuantos dejaron log  (ultimo job ID)
#   3. Rastros de error en los .err
#   4. Muestras con entrada vs muestras con salida
#   5. Salidas vacias o sospechosamente pequenas
#   + memoria y tiempo reales, para afinar el -M y el -W del siguiente paso
#
# Uso:
#   ./check_step.sh fastqc
#   ./check_step.sh cat
#   ./check_step.sh --all        # estado de todo el pipeline de un vistazo
#   ./check_step.sh --list
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

SAMPLE_FIELDS="${SAMPLE_FIELDS:-2}"     # M11_S3 -> 2 campos separados por '_'

STEPS_ALL=(cat fastqc polyg trim btmap unmap star spades quast genomad checkv diamond)

usage() {
  cat <<'EOF'
Uso: ./check_step.sh <paso> | --all | --list  [--exclude M41_S4,...]

Pasos:
  cat       Concatenado de lanes (compress.sh)
  fastqc    FastQC sobre cat_fastq
  polyg     Recorte de colas polyG (fastp)
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
}

STEPS=(); EXCLUDE=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --list|-h|--help) usage; exit 0 ;;
    --all)     STEPS=("${STEPS_ALL[@]}") ;;
    --exclude) EXCLUDE="${2:?--exclude necesita una lista separada por comas}"; shift ;;
    -*) echo "Opcion desconocida: $1"; usage; exit 1 ;;
    *)  STEPS+=("$1") ;;
  esac
  shift
done
[[ ${#STEPS[@]} -eq 0 ]] && { usage; exit 0; }

TOT_ERR=0; TOT_WARN=0

hum() { numfmt --to=iec-i --suffix=B "$1" 2>/dev/null || echo "${1}B"; }

# ---------------------------------------------------------------------------
# Configuracion por paso.
#
# in_key/out_key indican de donde sale el nombre de muestra:
#   file = del nombre del fichero    (M1_S1.bam -> M1_S1)
#   dir  = del nombre de su carpeta  (spades/M1_S1/contigs.fasta -> M1_S1)
#
# list_file = el .txt que se paso al job array. Sirve para saber cuantos
#             elementos DEBERIAN haber corrido, no solo cuantos corrieron.
# also_pat  = segundo patron de salida que tambien debe existir (p.ej. los R2).
# note      = aviso que se imprime antes de nada (dependencias no implementadas).
# ---------------------------------------------------------------------------
configure_step() {
  logpre=""; in_dir=""; in_pat=""; in_key=file
  out_dir=""; out_pat=""; out_key=file
  min_bytes=0; list_file=""; also_pat=""; note=""

  case "$1" in
    cat)
      logpre="compress"
      in_dir="$DATA_DIR";                      in_pat="*_R1_*.fastq.gz"
      out_dir="$SCRATCH/cat_fastq";            out_pat="*_R1_001.fastq.gz"
      also_pat="*_R2_001.fastq.gz"
      min_bytes=1000000 ;;
    fastqc)
      logpre="fastqc"
      in_dir="$SCRATCH/cat_fastq";             in_pat="*_R1_*.fastq.gz"
      out_dir="$RESULTS_DIR/fastqc_raw";       out_pat="*_R1*_fastqc.zip"
      also_pat="*_R2*_fastqc.zip"
      list_file="$SCRATCH/filelist.txt"
      min_bytes=100000 ;;
      polyg)
      logpre="polyg"
      in_dir="$SCRATCH/cat_fastq";             in_pat="*_R1_001.fastq.gz"
      out_dir="$SCRATCH/nopolyg";              out_pat="*_R1_001.fastq.gz"
      also_pat="*_R2_001.fastq.gz"
      list_file="$SCRATCH/samples.txt"
      min_bytes=1000000 ;;
    trim)
      logpre="trim"
      in_dir="$SCRATCH/cat_fastq";             in_pat="*_R1_*.fastq.gz" # in_dir="$SCRATCH/nopolyg";               in_pat="*_R1_*.fastq.gz"
      out_dir="$SCRATCH/trimmed/$RESULT_FROM"; out_pat="*_R1_001_paired.fastq.gz"
      also_pat="*_R2_001_paired.fastq.gz"
      list_file="$SCRATCH/samples.txt"
      min_bytes=10000000 ;;
    btmap)
      logpre="btmap"
      in_dir="$SCRATCH/trimmed/$RESULT_FROM";  in_pat="*_R1_001_paired.fastq.gz"
      out_dir="$SCRATCH/mapped/$RESULT_FROM";  out_pat="*.bam"
      list_file="$SCRATCH/samples.txt"
      min_bytes=10000000 ;;
    unmap)
      logpre="unmap"
      in_dir="$SCRATCH/mapped/$RESULT_FROM";   in_pat="*.bam"
      out_dir="$SCRATCH/unmapped_fastq";       out_pat="*_unmapped_R1.fastq.gz"
      also_pat="*_unmapped_R2.fastq.gz"
      list_file="$SCRATCH/bam_list.txt"
      min_bytes=1000 ;;
    star)
      logpre="star"
      in_dir="$SCRATCH/trimmed/$RESULT_FROM";  in_pat="*_R1_001_paired.fastq.gz"
      out_dir="$SCRATCH/aligned";              out_pat="*_Aligned.out.bam"
      min_bytes=1000000 ;;
    spades)
      logpre="spades"
      in_dir="$SCRATCH/unmapped_fastq";        in_pat="*_unmapped_R1.fastq.gz"
      out_dir="$SCRATCH/spades";               out_pat="contigs.fasta";  out_key=dir
      min_bytes=1000 ;;
    quast)
      logpre="quast"
      in_dir="$SCRATCH/spades";                in_pat="contigs.fasta";   in_key=dir
      out_dir="$RESULTS_DIR/quast";            out_pat="report.tsv";     out_key=dir
      list_file="$SCRATCH/assembly_list.txt"
      min_bytes=100 ;;
    genomad)
      logpre="genomad"
      in_dir="$SCRATCH/spades";                in_pat="*_final.fasta"
      out_dir="$RESULTS_DIR/genomad";          out_pat="*_summary";      out_key=dir
      min_bytes=1
      note="ningun script genera '*_final.fasta'. Falta el renombrado post-SPAdes." ;;
    checkv)
      logpre="checkv"
      in_dir="$SCRATCH/spades";                in_pat="*_final.fasta"
      out_dir="$RESULTS_DIR/checkv";           out_pat="quality_summary.tsv"; out_key=dir
      min_bytes=100
      note="ningun script genera '*_final.fasta'. Falta el renombrado post-SPAdes." ;;
    diamond)
      logpre="dmnd"
      in_dir="$SCRATCH/spades";                in_pat="contigs_1000bp.fasta"; in_key=dir
      out_dir="$RESULTS_DIR/diamond";          out_pat="*_diamond_extended.tsv"
      min_bytes=1
      note="ningun script genera 'contigs_1000bp.fasta'. Falta el filtro por longitud >=1000 bp." ;;
    *)
      echo "Paso desconocido: $1  (usa --list)"; return 1 ;;
  esac
  return 0
}

# Extrae el nombre de muestra de una ruta.
sample_key() {
  local path="$1" mode="$2" b
  if [[ "$mode" == "dir" ]]; then
    basename "$(dirname "$path")"
  else
    b=$(basename "$path")
    b="${b%%.*}"                            # quita desde el primer punto
    echo "$b" | cut -d_ -f1-"$SAMPLE_FIELDS"
  fi
}

keys_of() {   # keys_of <dir> <pat> <mode>
  find "$1" -name "$2" 2>/dev/null | while IFS= read -r p; do
    sample_key "$p" "$3"
  done | sort -u
}

# ---------------------------------------------------------------------------
run_step() {
  local step="$1"
  configure_step "$step" || return 1

  local ERRORS=0 WARNS=0
  err()  { echo "  [ERROR] $*"; ERRORS=$((ERRORS+1)); }
  warn() { echo "  [warn ] $*"; WARNS=$((WARNS+1)); }
  ok()   { echo "  [ok   ] $*"; }

  echo "==========================================================="
  echo "check_step.sh  paso: $step     $(date '+%Y-%m-%d %H:%M')"
  echo "==========================================================="
  [[ -n "$note" ]] && { echo; echo "  [AVISO] $note"; }
  echo

  # --- 1. Logs de LSF -------------------------------------------------------
  echo "== 1. Logs de LSF ($LOGS_DIR/$logpre.*)"
  shopt -s nullglob
  local logs=( "$LOGS_DIR/$logpre".*.log )
  shopt -u nullglob

  # Los logs se acumulan entre relanzamientos, cada uno con su %J. Si al
  # relanzar 3 indices fallidos contaramos todos los logs, los 3 fallos
  # antiguos seguirian apareciendo para siempre. Se toma UN log por indice:
  # el del intento mas reciente (mayor job ID).
  local logs_eff=()
  if [[ ${#logs[@]} -gt 0 ]]; then
    mapfile -t logs_eff < <(
      printf '%s\n' "${logs[@]}" | awk -F/ '{
        f = $NF; n = split(f, a, ".")
        idx = (n >= 4 ? a[3] : 0); job = a[2]
        print idx "\t" job "\t" $0
      }' | sort -k1,1n -k2,2n | awk -F'\t' '{ last[$1] = $3 } END { for (k in last) print last[k] }'
    )
  fi

  local n_ok=0 n_exit=0 n_mem=0 n_run=0
  if [[ ${#logs_eff[@]} -eq 0 ]]; then
    warn "no hay logs con prefijo '$logpre' (paso no lanzado todavia?)"
  else
    echo "  ${#logs[@]} logs en disco -> ${#logs_eff[@]} indices distintos (ultimo intento de cada uno)"
    n_ok=$(grep -l "Successfully completed"  "${logs_eff[@]}" 2>/dev/null | wc -l)
    n_exit=$(grep -l "Exited with"           "${logs_eff[@]}" 2>/dev/null | wc -l)
    n_mem=$(grep -l "TERM_MEMLIMIT"          "${logs_eff[@]}" 2>/dev/null | wc -l)
    n_run=$(grep -l "TERM_RUNLIMIT"          "${logs_eff[@]}" 2>/dev/null | wc -l)

    echo "  completados OK : $n_ok"
    [[ $n_exit -gt 0 ]] && err "salieron con error: $n_exit"
    [[ $n_mem  -gt 0 ]] && err "matados por MEMORIA (TERM_MEMLIMIT): $n_mem  -> sube el -M"
    [[ $n_run  -gt 0 ]] && err "matados por TIEMPO (TERM_RUNLIMIT): $n_run  -> sube el -W o cola mas larga"

    if [[ $n_exit -gt 0 || $n_mem -gt 0 || $n_run -gt 0 ]]; then
      echo "  elementos fallidos (indices a relanzar):"
      grep -l -e "Exited with" -e "TERM_MEMLIMIT" -e "TERM_RUNLIMIT" "${logs_eff[@]}" 2>/dev/null \
        | xargs -n1 basename 2>/dev/null | awk -F. '{print "         indice "$3"   (job "$2")"}' | sort -t' ' -k2n | head -10
    fi
  fi
  echo

  # --- 2. Elementos lanzados vs elementos con log ---------------------------
  echo "== 2. Cobertura del array"
  if [[ ${#logs_eff[@]} -eq 0 ]]; then
    echo "  (sin logs)"
  else
    local n_eff=${#logs_eff[@]} expected last_job
    last_job=$(printf '%s\n' "${logs[@]}" | xargs -n1 basename \
               | awk -F. '{print $2}' | sort -un | tail -1)
    echo "  indices con log: $n_eff   (ultimo job: $last_job)"

    if [[ -n "$list_file" && -f "$list_file" ]]; then
      expected=$(grep -cve '^[[:space:]]*$' "$list_file")
      echo "  esperados ($(basename "$list_file")): $expected"
      if   (( n_eff < expected )); then
        err "faltan $((expected - n_eff)) elementos sin log: no arrancaron o siguen en cola"
        echo "         comprueba con: bjobs -A $last_job"
      elif (( n_eff > expected )); then
        warn "hay mas indices con log que lineas en la lista (la lista cambio desde el envio?)"
      else
        ok "los $expected elementos dejaron log"
      fi
    else
      [[ -n "$list_file" ]] && warn "no existe $list_file: no se puede saber cuantos se lanzaron"
    fi
  fi
  echo

  # --- 3. Rastros de error en los .err --------------------------------------
  # Muchas herramientas escriben a stderr en condiciones normales (bowtie2 y
  # Trimmomatic vuelcan ahi sus estadisticas), asi que un .err con contenido
  # NO es un fallo. Se buscan firmas de error de verdad.
  echo "== 3. Rastros de error en los .err"
  local errpat='Traceback|Exception|command not found|No such file|Killed|Segmentation fault|CondaError|Could not find|error:|ERROR:'
  shopt -s nullglob
  local errs=( "$LOGS_DIR/$logpre".*.err )
  shopt -u nullglob
  if [[ ${#errs[@]} -eq 0 ]]; then
    echo "  (no hay ficheros .err)"
  else
    local hits
    hits=$(grep -lE "$errpat" "${errs[@]}" 2>/dev/null | wc -l)
    if (( hits > 0 )); then
      warn "$hits ficheros .err con rastros de error:"
      grep -lE "$errpat" "${errs[@]}" 2>/dev/null | sed 's/^/         /' | head -5
      echo "         revisalos: grep -nE '$errpat' <fichero>"
    else
      ok "${#errs[@]} ficheros .err sin firmas de error"
    fi
  fi
  echo

  # --- 4. Entradas vs salidas ----------------------------------------------
  echo "== 4. Muestras con entrada vs muestras con salida"
  local in_keys out_keys n_in n_out
  mapfile -t in_keys  < <(keys_of "$in_dir"  "$in_pat"  "$in_key")
  mapfile -t out_keys < <(keys_of "$out_dir" "$out_pat" "$out_key")
  # Muestras excluidas a proposito (p.ej. M41_S4, fallida en secuenciacion):
  # siguen teniendo entrada en disco pero no deben exigir salida.
  if [[ -n "$EXCLUDE" ]]; then
    local keep=() x skip
    for k in ${in_keys[@]+"${in_keys[@]}"}; do
      skip=0
      IFS=',' read -ra _exc <<< "$EXCLUDE"
      for x in "${_exc[@]}"; do [[ "$k" == "$x" ]] && skip=1; done
      (( skip )) || keep+=("$k")
    done
    in_keys=(${keep[@]+"${keep[@]}"})
  fi
  n_in=${#in_keys[@]}; n_out=${#out_keys[@]}

  echo "  entradas ($in_dir / $in_pat):  $n_in muestras"
  echo "  salidas  ($out_dir / $out_pat): $n_out muestras"

  if [[ $n_in -eq 0 ]]; then
    warn "no hay entradas: el paso anterior no ha corrido, o la ruta esta mal"
  else
    local missing extra
    missing=$(comm -23 <(printf '%s\n' "${in_keys[@]}") \
                       <(printf '%s\n' ${out_keys[@]+"${out_keys[@]}"}) )
    extra=$(comm -13   <(printf '%s\n' "${in_keys[@]}") \
                       <(printf '%s\n' ${out_keys[@]+"${out_keys[@]}"}) )
    if [[ -n "$missing" ]]; then
      # Sin ningun log, "falta la salida" no es un fallo: es que no se ha
      # lanzado. Solo es error si el paso corrio y aun asi falta gente.
      if (( ${#logs_eff[@]} == 0 )); then
        warn "$(printf '%s\n' "$missing" | wc -l) muestras pendientes (el paso no se ha lanzado)"
      else
        err "$(printf '%s\n' "$missing" | wc -l) muestras SIN salida:"
        printf '%s\n' "$missing" | sed 's/^/         /' | head -20
      fi
    else
      ok "todas las $n_in muestras tienen salida"
    fi
    if [[ -n "$extra" ]]; then
      warn "salidas sin entrada correspondiente (restos de otra corrida?):"
      printf '%s\n' "$extra" | sed 's/^/         /' | head -10
    fi
  fi

  # Segundo patron: los R2, que si no se cuentan pasan desapercibidos.
  if [[ -n "$also_pat" ]]; then
    local n_a n_p
    n_a=$(find "$out_dir" -name "$also_pat" 2>/dev/null | wc -l)
    n_p=$(find "$out_dir" -name "$out_pat"  2>/dev/null | wc -l)
    if (( n_a == n_p )); then
      ok "$out_pat: $n_p   |   $also_pat: $n_a   (cuadran)"
    else
      err "$out_pat: $n_p pero $also_pat: $n_a  -> falta media pareja"
    fi
  fi
  echo

  # --- 5. Salidas vacias o minusculas ---------------------------------------
  echo "== 5. Salidas vacias o sospechosas (< $(hum "$min_bytes"))"
  local small=() sz
  while IFS= read -r f; do
    [[ -z "$f" ]] && continue
    sz=$(stat -c %s "$f" 2>/dev/null || echo 0)
    (( sz < min_bytes )) && small+=("$(hum "$sz")  $(basename "$(dirname "$f")")/$(basename "$f")")
  done < <(find "$out_dir" \( -name "$out_pat" ${also_pat:+-o -name "$also_pat"} \) -type f 2>/dev/null)

  if (( ${#small[@]} == 0 )); then
    ok "ninguna salida vacia o minuscula"
  else
    warn "${#small[@]} salidas por debajo del umbral:"
    printf '         %s\n' "${small[@]}" | head -10
    (( ${#small[@]} > 10 )) && echo "         ... y $(( ${#small[@]} - 10 )) mas"
  fi
  echo

  # --- 6. Consumo real de memoria y tiempo ----------------------------------
  echo "== 6. Recursos usados realmente"
  if [[ ${#logs_eff[@]} -gt 0 ]]; then
    grep -h "Max Memory" "${logs_eff[@]}" 2>/dev/null \
      | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) {print $i; break}}' \
      | sort -n | awk '
          {v[NR]=$1}
          END{
            if(NR==0){print "  memoria: (LSF no reporto Max Memory)"; exit}
            printf "  memoria: min %s MB | mediana %s MB | max %s MB  (n=%d)\n", v[1], v[int((NR+1)/2)], v[NR], NR
            printf "  sugerencia -M: %d  (max x1.3)\n", v[NR]*1.3
          }'
    grep -h "Run time" "${logs_eff[@]}" 2>/dev/null \
      | awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/) {print $i; break}}' \
      | sort -n | awk '
          {v[NR]=$1}
          END{
            if(NR==0) exit
            printf "  tiempo : min %ds | mediana %ds | max %ds\n", v[1], v[int((NR+1)/2)], v[NR]
          }'
  else
    echo "  (sin logs)"
  fi
  echo

  # --- Resumen del paso -----------------------------------------------------
  echo "-----------------------------------------------------------"
  echo "paso '$step': $ERRORS errores, $WARNS avisos"
  if [[ ${#logs_eff[@]} -eq 0 && $n_out -eq 0 ]]; then
    echo "Este paso NO se ha ejecutado todavia."
  elif (( ERRORS > 0 )); then
    echo "Relanza SOLO los elementos que fallaron antes de continuar."
    echo "Truco: bsub -J \"${logpre}rerun[i,j,k]\" ...  (indices concretos entre corchetes)"
  elif (( WARNS > 0 )); then
    echo "Sin errores, pero revisa los avisos antes de continuar."
  else
    echo "Paso completo. Puedes continuar."
  fi
  echo "-----------------------------------------------------------"
  echo

  TOT_ERR=$((TOT_ERR + ERRORS)); TOT_WARN=$((TOT_WARN + WARNS))
  return 0
}

# ---------------------------------------------------------------------------
for s in "${STEPS[@]}"; do
  run_step "$s" || exit 1
done

if [[ ${#STEPS[@]} -gt 1 ]]; then
  echo "==========================================================="
  echo "TOTAL ${#STEPS[@]} pasos: $TOT_ERR errores, $TOT_WARN avisos"
  echo "==========================================================="
fi

exit $(( TOT_ERR > 0 ? 1 : 0 ))
