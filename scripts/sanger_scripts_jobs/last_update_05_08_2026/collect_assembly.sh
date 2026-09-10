#!/bin/bash
# ---------------------------------------------------------------------------
# collect_assembly.sh — reúne todo lo que hace falta para describir los
# ensamblajes: longitud y cobertura de cada contig, aporte de cada ensamblador
# según CD-HIT, y coste computacional.
#
# La longitud SIEMPRE se cuenta de la secuencia, no de la cabecera: SPAdes,
# MEGAHIT y Trinity las escriben de tres formas distintas y el `length` de la
# cabecera del union viene del ensamblador de origen, no del fichero.
#
# SALIDAS  ($RESULTS_DIR/qc_control/)
#   assembly_contigs.tsv.gz    set · sample · contig · origin · length · coverage
#   assembly_summary.tsv       set · sample · contigs · contigs_1kb · total_bp · longest · shortest
#   derep_clusters.tsv         sample · cluster · rep_len · n_members · assemblers
#   assembly_resources.tsv     tool · mode · sample · jobid · max_mem_mb · run_time_s
#   assembly_params.tsv        set · sample · tool · version · kmers · min_contig · extra · fuente
#
# USO
#   ./collect_assembly.sh --check     # qué encuentra y qué falta. EMPEZAR AQUI.
#   ./collect_assembly.sh             # escribe los cuatro ficheros
#
# Variables de entorno útiles:
#   MIN_LEN=500        longitud mínima de contig que se exporta (por defecto 500)
#   LISTS_DIR          carpeta con los fasta_<conjunto>.txt de make_fasta_lists.sh
#   DEREP_DIR          carpeta con los .clstr de CD-HIT
#   SAMPLES            lista de muestras
#   SET_ALIASES        renombrar conjuntos:  "spades=rnaviral"
#   SET_EXCLUDE        excluir conjuntos:    "union_all,trinity_rescate"
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

CHECK=0
[[ "${1:-}" == "--check" ]] && CHECK=1

MIN_LEN="${MIN_LEN:-500}"
OUT="$RESULTS_DIR/qc_control"; mkdir -p "$OUT"

CONTIGS="$OUT/assembly_contigs.tsv"
SUMMARY="$OUT/assembly_summary.tsv"
CLUSTERS="$OUT/derep_clusters.tsv"
RESOURCES="$OUT/assembly_resources.tsv"
PARAMS="$OUT/assembly_params.tsv"

first_dir() { local d; for d in "$@"; do [[ -d "$d" ]] && { echo "$d"; return; }; done; }

LISTS_DIR="${LISTS_DIR:-$(first_dir "$SCRATCH/lists")}"
DEREP_DIR="${DEREP_DIR:-$(first_dir "$SCRATCH/union" "$RESULTS_DIR/derep_stats" "$SCRATCH/derep")}"

# --- muestras ---------------------------------------------------------------
SAMPLES="${SAMPLES:-}"
if [[ -z "$SAMPLES" ]]; then
  for c in "$SCRATCH/asm_samples.txt" "$SCRATCH/samples.txt"; do
    [[ -s "$c" ]] && { SAMPLES="$c"; break; }
  done
fi
[[ -s "${SAMPLES:-}" ]] || { echo "ERROR: pasa SAMPLES=/ruta/lista.txt" >&2; exit 1; }
mapfile -t SAMPLE_LIST < <(grep -v '^[[:space:]]*$' "$SAMPLES" | sort -u)
N=${#SAMPLE_LIST[@]}

# --- conjuntos de ensamblaje ------------------------------------------------
# Fuente de verdad: los fasta_<conjunto>.txt que ya genera make_fasta_lists.sh.
# Si no están, se cae a un escaneo de directorios de $SCRATCH.
# Se juntan las dos fuentes: los fasta_<conjunto>.txt de make_fasta_lists.sh y
# un escaneo de directorios de $SCRATCH. Con solo las listas se pierden los
# modos que nunca llegaron a listarse.
declare -A SEEN=()
declare -A ALIAS_OF=()
declare -a SETS=()
add_set() { [[ -z "${SEEN[$1]:-}" ]] && { SEEN[$1]=1; SETS+=("$1"); }; }

if [[ -n "${LISTS_DIR:-}" && -d "$LISTS_DIR" ]]; then
  while IFS= read -r f; do
    [[ -s "$f" ]] || continue
    b=$(basename "$f"); b="${b#fasta_}"; b="${b%.txt}"
    add_set "$b"
  done < <(find "$LISTS_DIR" -maxdepth 1 -name 'fasta_*.txt' 2>/dev/null | sort)
fi

for d in "$SCRATCH"/spades* "$SCRATCH"/megahit* "$SCRATCH"/trinity* "$SCRATCH"/union*; do
  [[ -d "$d" ]] || continue
  b=$(basename "$d")
  # $SCRATCH/spades_rna -> "rna";  $SCRATCH/spades -> "spades" (rama rnaviral)
  b="${b#spades_}"
  add_set "$b"
done

# Exclusión opcional:  SET_EXCLUDE="union_all,trinity_rescate"
if [[ -n "${SET_EXCLUDE:-}" ]]; then
  IFS=',' read -ra _ex <<< "$SET_EXCLUDE"
  declare -a KEEP=()
  for x in "${SETS[@]}"; do
    skip=0
    for e in "${_ex[@]}"; do [[ "$x" == "$e" ]] && skip=1; done
    (( skip )) || KEEP+=("$x")
  done
  SETS=("${KEEP[@]}")
fi

# Renombrado opcional:  SET_ALIASES="spades=rnaviral,megahit_v2=megahit"
if [[ -n "${SET_ALIASES:-}" ]]; then
  IFS=',' read -ra _al <<< "$SET_ALIASES"
  for a in "${_al[@]}"; do
    from="${a%%=*}"; to="${a#*=}"
    for i in "${!SETS[@]}"; do [[ "${SETS[$i]}" == "$from" ]] && SETS[$i]="$to"; done
    ALIAS_OF["$to"]="$from"
  done
fi

# Directorio de un conjunto. Con lista, se deduce del primer fasta que contenga.
set_dir() {
  # OJO: bash expande TODOS los argumentos de `local` antes de asignar ninguno,
  # así que `local s="$1" lf="...${s}..."` deja $s sin definir y con set -u
  # aborta. Por eso van en líneas separadas.
  local s="$1"
  local lf="${LISTS_DIR:-/nonexistent}/fasta_${s}.txt"
  local first
  if [[ -s "$lf" ]]; then
    first=$(head -1 "$lf")
    [[ -n "$first" ]] && { dirname "$(dirname "$first")"; return; }
  fi
  local orig="${ALIAS_OF[$s]:-$s}"
  for d in "$SCRATCH/$s" "$SCRATCH/spades_$s" "$SCRATCH/$orig" "$SCRATCH/spades_$orig"; do
    [[ -d "$d" ]] && { echo "$d"; return; }
  done
}

# Fasta de una muestra dentro de un conjunto: SIEMPRE contigs.fasta, que es el
# ensamblaje completo. Los *_final.fasta de las listas ya vienen filtrados a
# >=1000 bp y con ellos no se puede dibujar la curva Nx ni el corte de 500.
sample_fasta() {
  local sdir="$1" s="$2" f
  for f in "$sdir/$s/contigs.fasta" "$sdir/$s/transcripts.fasta" "$sdir/$s/${s}_final.fasta"; do
    [[ -s "$f" ]] && { echo "$f"; return; }
  done
}

# --- parseo de un fasta -----------------------------------------------------
# La longitud se cuenta de la secuencia. De la cabecera solo se sacan la
# cobertura y, en el union, el ensamblador de origen (prefijo antes de "__").
#   SPAdes   >NODE_1_length_12898_cov_47.312
#   MEGAHIT  >k141_1234 flag=1 multi=47.3120 len=12898
#   Trinity  >TRINITY_DN1000_c0_g1_i1 len=329 path=[...]     (sin cobertura)
#   union    >metaviral__NODE_1_length_15437_cov_5.7
parse_fasta() {
  local f="$1" set="$2" s="$3"
  awk -v SET="$set" -v SAMPLE="$s" -v MINLEN="$MIN_LEN" '
    function flush(  cov, org, h, p) {
      if (name == "") return
      if (len >= MINLEN) {
        h = name; org = "NA"; cov = "NA"
        p = index(h, "__")
        if (p > 0) { org = substr(h, 1, p-1); h = substr(h, p+2) }
        if (match(h, /_cov_[0-9.]+/))       cov = substr(h, RSTART+5, RLENGTH-5)
        else if (match(h, /multi=[0-9.]+/)) cov = substr(h, RSTART+6, RLENGTH-6)
        sub(/[.]$/, "", cov)
        printf "%s\t%s\t%s\t%s\t%d\t%s\n", SET, SAMPLE, name, org, len, cov
      }
      name = ""; len = 0
    }
    /^>/ { flush(); name = substr($1, 2); len = 0; next }
    { len += length($0) }
    END { flush() }
  ' "$f"
}

# --- .clstr de CD-HIT -------------------------------------------------------
# Formato:
#   >Cluster 0
#   0	15437nt, >metaviral__NODE_1... *          <- el representante
#   1	15436nt, >megahit__k141_... at -/100.00%
parse_clstr() {
  local f="$1" s="$2"
  awk -v SAMPLE="$s" '
    function flush(  k, out) {
      if (cluster == "") return
      out = ""
      for (k in seen) out = (out == "" ? k : out "," k)
      printf "%s\t%s\t%d\t%s\t%d\t%s\n", SAMPLE, cluster, replen, repasm, nmem, (out == "" ? "NA" : out)
      delete seen
    }
    /^>Cluster/ { flush(); cluster = $2; replen = 0; repasm = "NA"; nmem = 0; next }
    {
      nmem++
      n = split($0, a, "nt, >")
      if (n < 2) next
      # a[1] acaba en la longitud; a[2] empieza por la cabecera
      m = split(a[1], b, "\t"); L = b[m] + 0
      hdr = a[2]
      p = index(hdr, "__")
      asm = (p > 0) ? substr(hdr, 1, p-1) : "unknown"
      seen[asm] = 1
      if ($0 ~ /\*[[:space:]]*$/) { replen = L; repasm = asm }
    }
    END { flush() }
  ' "$f"
}

# --- parámetros del ensamblador --------------------------------------------
# De los ficheros que el propio ensamblador deja en su directorio de salida.
# NADA se da por sabido: si no aparece en un fichero, sale "NA" y se ve.
#   SPAdes   spades.log / params.txt   -> "k values to be used: 33, 49"
#   MEGAHIT  log / opts.txt            -> "k list: 21,29,39,..."
#   Trinity  Trinity.timing / el log   -> k fijo a 25, no configurable
parse_params() {
  local sdir="$1" set="$2" s="$3"
  local d="$sdir/$s" tool="NA" ver="NA" kms="NA" minc="NA" extra="NA" src="NA"
  local f

  # ---- SPAdes
  for f in "$d/spades.log" "$d/params.txt"; do
    [[ -s "$f" ]] || continue
    tool="SPAdes"; src="$(basename "$f")"
    [[ "$ver" == "NA" ]] && ver=$(grep -om1 'SPAdes v[0-9.]*' "$f" 2>/dev/null | head -1)
    [[ -z "$ver" ]] && ver="NA"
    if [[ "$kms" == "NA" ]]; then
      kms=$(grep -oiE 'k *values? *(to be used)? *:? *\[?[0-9]+(, *[0-9]+)*' "$f" 2>/dev/null \
            | head -1 | grep -oE '[0-9]+(, *[0-9]+)*' | tr -d ' ')
      [[ -z "$kms" ]] && kms=$(grep -oE 'with K=[0-9]+' "$f" 2>/dev/null \
            | grep -oE '[0-9]+' | sort -n -u | paste -sd, -)
      [[ -z "$kms" ]] && kms="NA"
    fi
    [[ "$extra" == "NA" ]] && extra=$(grep -om1 -- '--\(rnaviral\|metaviral\|meta\|rna\|isolate\|careful\)' "$f" 2>/dev/null | tr -d '-')
    [[ -z "$extra" ]] && extra="NA"
  done

  # ---- MEGAHIT
  if [[ "$tool" == "NA" ]]; then
    for f in "$d/log" "$d/opts.txt" "$d"/*.log; do
      [[ -s "$f" ]] || continue
      grep -qi megahit "$f" 2>/dev/null || continue
      tool="MEGAHIT"; src="$(basename "$f")"
      ver=$(grep -om1 'MEGAHIT v[0-9.]*' "$f" 2>/dev/null | head -1); [[ -z "$ver" ]] && ver="NA"
      kms=$(grep -oiE 'k list: *[0-9,]+' "$f" 2>/dev/null | head -1 | grep -oE '[0-9,]+$')
      [[ -z "$kms" ]] && kms="NA"
      minc=$(grep -oE 'min_contig_len[ =:]+[0-9]+' "$f" 2>/dev/null | head -1 | grep -oE '[0-9]+$')
      [[ -z "$minc" ]] && minc="NA"
      break
    done
  fi

  # ---- Trinity
  if [[ "$tool" == "NA" ]] && compgen -G "$d/Trinity*" >/dev/null 2>&1; then
    tool="Trinity"; src="directorio Trinity"
    for f in "$d"/Trinity.timing "$d"/*.log; do
      [[ -s "$f" ]] || continue
      ver=$(grep -om1 'Trinity-v[0-9.]*' "$f" 2>/dev/null | head -1); [[ -n "$ver" ]] && break
    done
    [[ -z "$ver" ]] && ver="NA"
    kms="25"                # Trinity NO permite cambiarlo
    minc=$(grep -oE 'min_contig_length[ =]+[0-9]+' "$d"/*.log 2>/dev/null | head -1 | grep -oE '[0-9]+$')
    [[ -z "$minc" ]] && minc="NA"
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$set" "$s" "$tool" "$ver" "$kms" "$minc" "$extra" "$src"
}

# --- recursos de LSF --------------------------------------------------------
# "Max Memory :  968 MB" y "Run time :  852 sec." del pie de los logs de LSF.
parse_lsf() {
  local f="$1"
  awk '
    /Max Memory[[:space:]]*:/ { for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) { mem=$i; unit=$(i+1); break } }
    /Run time[[:space:]]*:/   { for (i=1;i<=NF;i++) if ($i ~ /^[0-9]+$/) { rt=$i; break } }
    END {
      if (unit == "GB") mem = mem * 1024
      if (mem == "" ) mem = "NA"
      if (rt  == "" ) rt  = "NA"
      printf "%s\t%s\n", mem, rt
    }
  ' "$f"
}

# ---------------------------------------------------------------------------
# --check
# ---------------------------------------------------------------------------
if (( CHECK )); then
  printf '\n== Fuentes ==\n'
  printf '  %-20s %s\n' "muestras"  "$SAMPLES ($N)"
  printf '  %-20s %s\n' "listas"    "${LISTS_DIR:-[NO ENCONTRADO]}"
  printf '  %-20s %s\n' "derep"     "${DEREP_DIR:-[NO ENCONTRADO]}"
  printf '  %-20s %s\n' "logs LSF"  "$LOGS_DIR"
  printf '  %-20s %s\n' "MIN_LEN"   "$MIN_LEN bp"

  printf '\n== Conjuntos de ensamblaje (%d) ==\n' "${#SETS[@]}"
  if [[ ${#SETS[@]} -eq 0 ]]; then
    printf '  ninguno. Corre make_fasta_lists.sh, o pasa LISTS_DIR=...\n'
  fi
  for set in "${SETS[@]}"; do
    sdir=$(set_dir "$set"); hit=0
    if [[ -n "$sdir" ]]; then
      for s in "${SAMPLE_LIST[@]}"; do
        [[ -n "$(sample_fasta "$sdir" "$s")" ]] && hit=$((hit+1))
      done
    fi
    printf '  %-14s %2d/%d muestras   %s\n' "$set" "$hit" "$N" "${sdir:-[sin directorio]}"
  done

  cl=0
  for s in "${SAMPLE_LIST[@]}"; do
    [[ -n "$(find "${DEREP_DIR:-/nonexistent}" -maxdepth 2 -name "${s}*.clstr" 2>/dev/null | head -1)" ]] && cl=$((cl+1))
  done
  printf '\n== CD-HIT ==\n  ficheros .clstr encontrados: %d/%d\n' "$cl" "$N"
  if (( cl == 0 )); then
    printf '  Busca dónde están y pásalo:  find "$SCRATCH" -name "*.clstr" | head -3\n'
    find "$SCRATCH" -maxdepth 4 -name "*.clstr" 2>/dev/null | head -3 | sed 's/^/    /'
  fi

  nlog=$(find "$LOGS_DIR" -maxdepth 1 \( -name 'spades*.log' -o -name 'megahit*.log' -o -name 'trinity*.log' \) 2>/dev/null | wc -l)
  printf '\n== Recursos ==\n  logs de ensamblaje: %d\n' "$nlog"
  (( nlog == 0 )) && printf '  Sin logs no hay figura de coste; el resto funciona igual.\n'
  printf '\n'
  exit 0
fi

# ---------------------------------------------------------------------------
# RECOLECCIÓN
# ---------------------------------------------------------------------------
printf 'set\tsample\tcontig\torigin\tlength\tcoverage\n' > "$CONTIGS"
printf 'set\tsample\tcontigs\tcontigs_1kb\ttotal_bp\tlongest\tshortest\n' > "$SUMMARY"
printf 'sample\tcluster\trep_len\trep_asm\tn_members\tassemblers\n' > "$CLUSTERS"
printf 'tool\tmode\tsample\tjobid\tmax_mem_mb\trun_time_s\n'    > "$RESOURCES"
printf 'set\tsample\ttool\tversion\tkmers\tmin_contig\textra\tfuente\n' > "$PARAMS"

echo "Conjuntos: ${SETS[*]}"
for set in "${SETS[@]}"; do
  sdir=$(set_dir "$set")
  [[ -n "$sdir" ]] || { echo "  $set: sin directorio, se salta"; continue; }
  n_ok=0
  for s in "${SAMPLE_LIST[@]}"; do
    f=$(sample_fasta "$sdir" "$s"); [[ -n "$f" ]] || continue
    parse_fasta "$f" "$set" "$s" >> "$CONTIGS"
    parse_params "$sdir" "$set" "$s" >> "$PARAMS"
    # Resumen sobre el fichero completo, sin filtro de MIN_LEN.
    awk -v SET="$set" -v SAMPLE="$s" '
      function acc() { if (n>0) { tot+=len; if (len>=1000) k++
                                  if (len>max) max=len
                                  if (min==0 || len<min) min=len } }
      /^>/ { acc(); n++; len=0; next }
      { len += length($0) }
      END { acc(); printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\n", SET, SAMPLE, n, k+0, tot+0, max+0, min+0 }
    ' "$f" >> "$SUMMARY"
    n_ok=$((n_ok+1))
  done
  echo "  $set: $n_ok/$N muestras"
done

echo "CD-HIT:"
n_cl=0
for s in "${SAMPLE_LIST[@]}"; do
  f=$(find "${DEREP_DIR:-/nonexistent}" -maxdepth 2 -name "${s}*.clstr" 2>/dev/null \
      | xargs -r ls -t 2>/dev/null | head -1)
  [[ -n "$f" ]] || continue
  parse_clstr "$f" "$s" >> "$CLUSTERS"
  n_cl=$((n_cl+1))
done
echo "  $n_cl/$N muestras con .clstr"

echo "Recursos de LSF:"
n_res=0
while IFS= read -r log; do
  [[ -s "$log" ]] || continue
  base=$(basename "$log")
  tool="${base%%.*}"
  jobid=$(echo "$base" | awk -F. '{print $2"."$3}')
  # modo: lo que diga el propio log
  mode=$(grep -om1 -- '--\(rnaviral\|metaviral\|meta\|rna\)\b' "$log" 2>/dev/null | tr -d '-')
  [[ -n "$mode" ]] || mode="NA"
  # muestra: la primera de la lista que aparezca en el log
  smp=$(grep -oFf <(printf '%s\n' "${SAMPLE_LIST[@]}") "$log" 2>/dev/null | head -1)
  [[ -n "$smp" ]] || smp="NA"
  read -r mem rt < <(parse_lsf "$log")
  [[ "$mem" == "NA" && "$rt" == "NA" ]] && continue
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$tool" "$mode" "$smp" "$jobid" "$mem" "$rt" >> "$RESOURCES"
  n_res=$((n_res+1))
done < <(find "$LOGS_DIR" -maxdepth 1 \
           \( -name 'spades*.log' -o -name 'megahit*.log' -o -name 'trinity*.log' \) 2>/dev/null | sort)
echo "  $n_res logs con recursos"

gzip -f "$CONTIGS"
echo
echo "Escrito en $OUT:"
ls -lh "$CONTIGS.gz" "$SUMMARY" "$CLUSTERS" "$RESOURCES" "$PARAMS" 2>/dev/null | awk '{printf "  %-10s %s\n", $5, $9}'
echo
echo "Cópialos a tu equipo:"
echo "  scp <usuario>@<login>:$OUT/{assembly_contigs.tsv.gz,assembly_summary.tsv,derep_clusters.tsv,assembly_resources.tsv,assembly_params.tsv} ."
