#!/bin/bash
# cribado_genoma.sh — que fraccion de cada vOTU esta en el genoma de referencia del huesped
#
# Trocea cada secuencia en ventanas de 100 pb (paso 50), las alinea con bowtie2 --local
# contra BT2_INDEX (AaegL5) y cuenta cuantas caen en el genoma casi identicas
# (>= 80 pb alineados y NM <= 2, es decir >= 98% de identidad).
#
# Una vOTU con frac_en_genoma >= 0.5 es secuencia del huesped: EVE/NIRVS transcrito
# si tiene hit viral, transcrito del mosquito si no. Detecta solo lo que esta en la
# referencia; un EVE polimorfico ausente de AaegL5 no aparece (eso lo dice el small RNA).
#
# Uso (LSF; el indice de Aedes necesita ~8 GB):
#   bsub -J cribado -q normal -n 4 -M 8000 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        -o "$LOGS_DIR/cribado.%J.log" -e "$LOGS_DIR/cribado.%J.err" \
#        "$SCRIPTS_DIR/cribado_genoma.sh"
#
# Salida: $RESULTS_DIR/qc_control/votus_en_genoma.tsv
#   votu  longitud  n_ventanas  n_alineadas  n_en_genoma  frac_en_genoma  cromosomas

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

IN="${1:-$RESULTS_DIR/qc_control/votus}"
OUT="${2:-$RESULTS_DIR/qc_control/votus_en_genoma.tsv}"
W=100; STEP=50; MAXNM=2; MINALN=80

[[ -s "$IN" ]] || { echo "ERROR: no existe $IN"; exit 1; }
[[ -n "${BT2_INDEX:-}" ]] || { echo "ERROR: BT2_INDEX sin definir"; exit 1; }

activate_env "${ENV_SAMTOOLS:-bowtie2_2.5.5}"
THREADS="${LSB_DJOB_NUMPROC:-4}"
export LC_ALL=C

tmp="$SCRATCH/tmp_cribado_${LSB_JOBID:-manual}"
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

echo "[cribado] entrada: $IN"
echo "[cribado] indice:  $BT2_INDEX"

# 1. ventanas (las secuencias de menos de W van enteras)
awk -v W="$W" -v S="$STEP" -v NV="$tmp/nvent.tsv" '
  function emit(   L, n, i) {
    if (id == "") return
    L = length(s); n = 0
    if (L <= W) { printf ">%s|1\n%s\n", id, s; n = 1 }
    else for (i = 1; i + W - 1 <= L; i += S) { printf ">%s|%d\n%s\n", id, i, substr(s, i, W); n++ }
    print id "\t" L "\t" n > NV
  }
  /^>/ { emit(); id = substr($1, 2); s = ""; next }
  { s = s $0 }
  END { emit() }' "$IN" > "$tmp/win.fa"
echo "[cribado] secuencias: $(wc -l < "$tmp/nvent.tsv")   ventanas: $(grep -c '^>' "$tmp/win.fa")"

# 2. alineamiento (mejor hit por ventana)
bowtie2 -p "$THREADS" -f --local -x "$BT2_INDEX" -U "$tmp/win.fa" \
        --no-unal --no-hd 2> "$tmp/bt2.log" > "$tmp/win.sam"
tail -1 "$tmp/bt2.log" | sed 's/^/[cribado] bowtie2: /'

# 3. resumen por vOTU
awk -F'\t' -v OFS='\t' -v W="$W" -v MAXNM="$MAXNM" -v MINALN="$MINALN" '
  NR == FNR { L[$1] = $2; n[$1] = $3; orden[++k] = $1; next }
  {
    v = $1; sub(/\|[0-9]+$/, "", v)
    cig = $6; aln = 0
    while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
      x = substr(cig, RSTART, RLENGTH); op = substr(x, length(x))
      if (op ~ /[M=X]/) aln += substr(x, 1, length(x) - 1) + 0
      cig = substr(cig, RSTART + RLENGTH)
    }
    nm = 99
    for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) nm = substr($i, 6) + 0
    minaln = (L[v] < W) ? int(0.8 * L[v]) : MINALN
    al[v]++
    if (aln >= minaln && nm <= MAXNM) {
      ok[v]++
      if (index("," chr[v] ",", "," $3 ",") == 0) chr[v] = (chr[v] == "" ? $3 : chr[v] "," $3)
    }
  }
  END {
    print "votu", "longitud", "n_ventanas", "n_alineadas", "n_en_genoma", "frac_en_genoma", "cromosomas"
    for (j = 1; j <= k; j++) {
      v = orden[j]
      printf "%s\t%d\t%d\t%d\t%d\t%.3f\t%s\n", v, L[v], n[v], al[v] + 0, ok[v] + 0, (ok[v] + 0) / n[v], chr[v]
    }
  }' "$tmp/nvent.tsv" "$tmp/win.sam" > "$OUT.partial"
mv "$OUT.partial" "$OUT"

# 4. informe
echo
echo "[cribado] distribucion de frac_en_genoma:"
awk -F'\t' 'NR > 1 { f = $6
  if (f == 0) b = "0"; else if (f < 0.25) b = "(0, 0.25)"; else if (f < 0.5) b = "[0.25, 0.5)"
  else if (f < 0.9) b = "[0.5, 0.9)"; else b = ">= 0.9"
  c[b]++ } END { split("0|(0, 0.25)|[0.25, 0.5)|[0.5, 0.9)|>= 0.9", o, "|")
  for (i = 1; i <= 5; i++) printf "    %-12s %4d\n", o[i], c[o[i]] }' "$OUT"
echo "[cribado] escrito: $OUT"
