#!/usr/bin/env bash
# =============================================================================
# isv_perfiles.sh -- Perfiles de cobertura por posicion sobre las referencias
#                    del chequeo de sanidad (isv_sanity.sh)
#
#   Un virus que se replica da cobertura pareja de extremo a extremo.
#   Un EVE o un repetido del hospedero da un parche corto sobre fondo cero.
#   Esta es la diferencia que ninguna tabla transmite y una figura si.
#
#   Lee los BAM que ya dejo isv_sanity.sh y emite una tabla larga, lista
#   para ggplot, con la profundidad media por ventana.
#
# Subcomandos:
#   check   verifica BAMs, entorno y referencias
#   run     calcula los perfiles de las 22 muestras (segundos; no necesita cola)
#   peek    imprime el perfil de una muestra/referencia para inspeccion rapida
#
# Variables utiles:
#   WIN=25        tamano de ventana en bp
#   BAM_TAG=.id   BAM filtrado por identidad (vacio = BAM crudo)
#   BASE_REF=     si se define, ademas emite resolucion de 1 bp para esa referencia
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${SCRIPT_DIR}/config.sh"

ISV_DIR="${REFS_DIR}/isv_sanity"
OUT_DIR="${RESULTS_DIR}/isv_sanity"
BAM_DIR="${OUT_DIR}/bam"
PERF_DIR="${OUT_DIR}/perfiles"
mkdir -p "${PERF_DIR}"

READ_SET="${READ_SET:-unmapped}"
BAM_TAG="${BAM_TAG-.id}"
WIN="${WIN:-25}"
BASE_REF="${BASE_REF:-}"
SAMPLES_FILE="${SAMPLES_FILE:-${ISV_DIR}/samples.txt}"

RED=$'\033[31m'; GRN=$'\033[32m'; YLW=$'\033[33m'; RST=$'\033[0m'
log()  { printf '%s[%s]%s %s\n' "${GRN}" "$(date +%H:%M:%S)" "${RST}" "$*"; }
warn() { printf '%s[WARN]%s %s\n' "${YLW}" "${RST}" "$*" >&2; }
die()  { printf '%s[ERROR]%s %s\n' "${RED}" "${RST}" "$*" >&2; exit 1; }

bam_de() { printf '%s/%s.%s%s.bam\n' "${BAM_DIR}" "$1" "${READ_SET}" "${BAM_TAG}"; }

# ---------------------------------------------------------------------------
cmd_check() {
  local fail=0
  log "Entorno ${ENV_BOWTIE2}:"
  if ( activate_env "${ENV_BOWTIE2}" >/dev/null 2>&1; command -v samtools >/dev/null ); then
    ( activate_env "${ENV_BOWTIE2}" >/dev/null 2>&1
      printf '  samtools   %s\n' "$(samtools --version | head -1 | awk '{print $NF}')" )
  else
    printf '  %ssamtools no disponible%s\n' "${RED}" "${RST}"; fail=1
  fi

  [[ -s "${SAMPLES_FILE}" ]] || die "falta ${SAMPLES_FILE}"
  log "BAMs (READ_SET=${READ_SET}, BAM_TAG='${BAM_TAG}'):"
  local n=0 falta=0
  while read -r s; do
    [[ -z "$s" ]] && continue
    n=$((n+1))
    local b; b="$(bam_de "$s")"
    if [[ -s "$b" && -s "${b}.bai" ]]; then :; else
      printf '  %sFALTA%s %s\n' "${RED}" "${RST}" "$b"; falta=$((falta+1))
    fi
  done < "${SAMPLES_FILE}"
  printf '  %d de %d presentes e indexados\n' "$((n-falta))" "$n"
  [[ $falta -eq 0 ]] || fail=1

  log "Referencias en el indice:"
  if [[ -s "${ISV_DIR}/isv_refs.fasta.fai" ]]; then
    awk -F'\t' '{printf "  %-30s %7d bp\n", $1, $2}' "${ISV_DIR}/isv_refs.fasta.fai"
  else
    printf '  %sfalta isv_refs.fasta.fai%s\n' "${RED}" "${RST}"; fail=1
  fi

  printf '  ventana = %s bp\n' "${WIN}"
  [[ $fail -eq 0 ]] && log "check OK" || die "check con fallos"
}

# ---------------------------------------------------------------------------
cmd_run() {
  [[ -s "${SAMPLES_FILE}" ]] || die "falta ${SAMPLES_FILE}"
  local out="${PERF_DIR}/perfiles_w${WIN}.${READ_SET}.tsv"
  local base_out="${PERF_DIR}/perfil_base_${BASE_REF}.${READ_SET}.tsv"

  printf 'muestra\tref\tvirus\tpos_ini\tpos_fin\tprof_media\tprof_max\tbases_cero\n' > "${out}"
  [[ -n "${BASE_REF}" ]] && printf 'muestra\tref\tpos\tprof\n' > "${base_out}"

  activate_env "${ENV_BOWTIE2}" >/dev/null

  local i=0
  while read -r s; do
    [[ -z "$s" ]] && continue
    i=$((i+1))
    local b; b="$(bam_de "$s")"
    [[ -s "$b" ]] || { warn "sin BAM: $s"; continue; }

    # -a  todas las posiciones, incluidas las de profundidad 0
    # -J  cuenta las deleciones dentro del alineamiento
    # -d 0 sin tope de profundidad (hay ventanas con decenas de miles de lecturas)
    samtools depth -a -J -d 0 "$b" \
      | awk -v S="$s" -v W="${WIN}" '
          { w = int(($2 - 1) / W) * W + 1
            k = $1 SUBSEP w
            if (!(k in n)) { ref[k] = $1; ini[k] = w }
            sum[k] += $3; n[k]++
            if ($3 > mx[k]) mx[k] = $3
            if ($3 == 0) z[k]++
          }
          END { for (k in n) {
                  split(ref[k], p, "__")
                  printf "%s\t%s\t%s\t%d\t%d\t%.3f\t%d\t%d\n",
                         S, ref[k], p[1], ini[k], ini[k] + n[k] - 1,
                         sum[k] / n[k], mx[k] + 0, z[k] + 0
                } }' \
      | sort -k2,2 -k4,4n >> "${out}"

    if [[ -n "${BASE_REF}" ]]; then
      samtools depth -a -J -d 0 -r "${BASE_REF}" "$b" \
        | awk -v S="$s" '{printf "%s\t%s\t%s\t%s\n", S, $1, $2, $3}' >> "${base_out}"
    fi
    printf '  %-10s hecho\n' "$s"
  done < "${SAMPLES_FILE}"

  log "${i} muestras -> ${out}  ($(wc -l < "${out}") filas)"
  [[ -n "${BASE_REF}" ]] && log "resolución de 1 bp -> ${base_out}"

  # Tamaño de librería, para normalizar entre muestras si hace falta
  local lib="${PERF_DIR}/tamanos_libreria.${READ_SET}.tsv"
  printf 'muestra\tpares_totales\tpares_alineados\ttasa\n' > "${lib}"
  while read -r s; do
    [[ -z "$s" ]] && continue
    local lg="${LOGS_DIR}/bt2_${s}.${READ_SET}.log"
    [[ -s "$lg" ]] || continue
    awk -v S="$s" '
      /reads; of these/        { tot = $1 }
      /overall alignment rate/ { gsub(/%/, "", $1); rate = $1 }
      END { printf "%s\t%d\t%d\t%.4f\n", S, tot, tot * rate / 100, rate }' "$lg" >> "${lib}"
  done < "${SAMPLES_FILE}"
  log "tamaños de librería -> ${lib}"

  log "Resumen: ventanas con prof_media ≥1 por virus"
  awk -F'\t' 'NR>1 && $6>=1 {c[$3]++} END{for(v in c) printf "  %-16s %6d ventanas\n", v, c[v]}' \
    "${out}" | sort
}

# ---------------------------------------------------------------------------
# peek -- perfil en texto, para ver la forma sin salir de la terminal
cmd_peek() {
  local s="${1:?uso: $0 peek <muestra> <ref>}"
  local r="${2:?uso: $0 peek <muestra> <ref>}"
  local out="${PERF_DIR}/perfiles_w${WIN}.${READ_SET}.tsv"
  [[ -s "${out}" ]] || die "no existe ${out}; corre '$0 run'"

  awk -F'\t' -v S="$s" -v R="$r" '
    $1==S && $2 ~ R { d[n] = $6; p[n] = $4; n++; if ($6 > mx) mx = $6 }
    END {
      if (n == 0) { print "sin datos para " S " / " R; exit 1 }
      printf "%s  %s   max %.0f×   %d ventanas\n\n", S, R, mx, n
      for (i = 0; i < n; i++) {
        h = (mx > 0) ? int(50 * d[i] / mx) : 0
        bar = ""
        for (j = 0; j < h; j++) bar = bar "#"
        printf "%7d |%-50s %8.1f\n", p[i], bar, d[i]
      }
    }' "${out}"
}

# ---------------------------------------------------------------------------
case "${1:-}" in
  check) cmd_check ;;
  run)   cmd_run ;;
  peek)  shift; cmd_peek "$@" ;;
  *) sed -n '2,24p' "${BASH_SOURCE[0]}"; exit 1 ;;
esac
