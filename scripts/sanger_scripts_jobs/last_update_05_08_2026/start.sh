#!/bin/bash
# ---------------------------------------------------------------------------
# start.sh — ritual de arranque del pipeline jkvirome (farm de Sanger)
#
# HAY QUE SOURCEARLO, NO EJECUTARLO:
#
#     source ./start.sh
#
# Si lo ejecutas (./start.sh) corre en un subshell y las variables mueren con
# el: no verias nada cargado. El script lo detecta y te avisa.
#
# Que hace:
#   1. Sourcea config.sh              (rutas y entornos)
#   2. Exporta LSB_DEFAULT_USERGROUP  (sin esto bsub aborta)
#   3. Comprueba que las rutas existen
#   4. Lista los entornos conda que existen de verdad
#   5. Cuota de Lustre y jobs en curso
#
# No hace ningun cambio: solo carga y reporta.
# ---------------------------------------------------------------------------

# --- 0. Debe sourcearse -----------------------------------------------------
if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo "ERROR: este script hay que sourcearlo, no ejecutarlo."
    echo "       source $(basename "${BASH_SOURCE[0]}")"
    exit 1
fi

# OJO: nada de 'set -e' aqui. Al sourcear, un fallo cerraria TU terminal.

_JKV_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

_ok()   { printf '  [ok  ] %s\n' "$*"; }
_warn() { printf '  [!!  ] %s\n' "$*"; }

echo "==========================================================="
echo "jkvirome — arranque    $(date '+%Y-%m-%d %H:%M')    $(hostname)"
echo "==========================================================="
echo

# --- 1. config.sh -----------------------------------------------------------
if [[ -f "$_JKV_DIR/config.sh" ]]; then
    source "$_JKV_DIR/config.sh"
    _ok "config.sh cargado desde $_JKV_DIR"
else
    _warn "NO existe $_JKV_DIR/config.sh — no puedo continuar"
    return 1
fi

# --- 2. Grupo LSF -----------------------------------------------------------
export LSB_DEFAULT_USERGROUP="${LSB_DEFAULT_USERGROUP:-team222}"
_ok "LSB_DEFAULT_USERGROUP=$LSB_DEFAULT_USERGROUP"

# Atajo para las lineas de bsub del RUNBOOK
S="$_JKV_DIR"
export S
_ok "\$S=$S   (usalo en bsub: \"\$S/fastqc_worker.sh ...\")"
echo

# --- 3. Rutas ---------------------------------------------------------------
echo "-- Rutas"
_missing=0
for v in DATA_DIR SCRATCH RESULTS_DIR LOGS_DIR REFS_DIR; do
    if [[ -d "${!v}" ]]; then
        printf '  [ok  ] %-12s %s\n' "$v" "${!v}"
    else
        printf '  [!!  ] %-12s %s   NO EXISTE\n' "$v" "${!v}"
        _missing=$((_missing + 1))
    fi
done
[[ ! -d "$DATA_DIR" ]] && _warn "DATA_DIR es el unico que no crea nadie: revisalo antes de seguir"
[[ ! -d "$REFS_DIR" ]] && _warn "REFS_DIR falta: bloquea trimming (ADAPTERS) y mapeo (BT2_INDEX)"
echo

# --- 4. Datos ---------------------------------------------------------------
echo "-- Datos en DATA_DIR"
if [[ -d "$DATA_DIR" ]]; then
    _n=$(find "$DATA_DIR" -type f -name "*.fastq.gz" 2>/dev/null | wc -l)
    _s=$(find "$DATA_DIR" -type f -name "*_R1_*.fastq.gz" -printf '%f\n' 2>/dev/null \
         | cut -d_ -f1-2 | sort -u | wc -l)
    _lanes=$(find "$DATA_DIR" -type f -name "*_L00*_R1_*.fastq.gz" 2>/dev/null | wc -l)
    printf '  %s ficheros .fastq.gz   |   %s muestras\n' "$_n" "$_s"
    if [[ $_lanes -gt 0 ]]; then
        _warn "hay lanes sin concatenar ($_lanes ficheros con _L00X_)"
        _warn "los workers esperan \${sample}_R1_001.fastq.gz -> pasa compress.sh antes del trimming"
    fi
else
    _warn "sin DATA_DIR no puedo contar"
fi
echo

# --- 5. Entornos conda ------------------------------------------------------
echo "-- Entornos conda"
if command -v conda >/dev/null 2>&1 || module load conda 2>/dev/null; then
    _envs=$(conda env list 2>/dev/null | awk '!/^#/ && NF {print $1}')
    for pair in "FASTQC:$ENV_FASTQC" "MULTIQC:$ENV_MULTIQC" "TRIMMOMATIC:$ENV_TRIMMOMATIC" \
                "BOWTIE2:$ENV_BOWTIE2" "SPADES:$ENV_SPADES" "QUAST:$ENV_QUAST" \
                "STAR:$ENV_STAR" "CHECKV:$ENV_CHECKV" "GENOMAD:$ENV_GENOMAD" \
                "DIAMOND:$ENV_DIAMOND" "COBRA:$ENV_COBRA"; do
        _label="${pair%%:*}"; _env="${pair#*:}"
        if grep -qx -- "$_env" <<< "$_envs"; then
            printf '  [ok  ] %-12s %s\n' "$_label" "$_env"
        else
            printf '  [!!  ] %-12s %s   NO EXISTE\n' "$_label" "$_env"
        fi
    done
else
    _warn "conda no disponible aqui: 'module load conda' para comprobar"
fi
echo

# --- 6. Cuota y jobs --------------------------------------------------------
echo "-- Cuota de Lustre"
lfs quota -h "$SCRATCH" 2>/dev/null | grep -E '^\s+[0-9]' | sed 's/^/  /' \
    || echo "  (lfs quota no disponible)"
echo

echo "-- Jobs en curso"
if command -v bjobs >/dev/null 2>&1; then
    bjobs -w 2>/dev/null | sed 's/^/  /' || echo "  ninguno"
else
    echo "  (bjobs no disponible: no estas en un nodo de login?)"
fi
echo

echo "==========================================================="
echo "Listo. Atajos:  \$S (scripts)  \$DATA_DIR  \$SCRATCH  \$RESULTS_DIR  \$LOGS_DIR"
echo "Control:        ./check_inputs.sh | ./check_step.sh <paso> | ./cleanup.sh"
echo "==========================================================="

unset _JKV_DIR _missing _n _s _lanes _envs _label _env pair
unset -f _ok _warn
