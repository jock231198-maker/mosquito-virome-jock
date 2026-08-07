#!/usr/bin/env bash
#
# recreate_all_envs.sh  --  EJECUTAR EN EL FARM (farm22)
#
# Recrea todos los entornos a partir de los YAML exportados desde el Mac.
# Los entornos que fallen se registran para revisarlos a mano.
#
# Uso:
#   chmod +x recreate_all_envs.sh
#   ./recreate_all_envs.sh ~/conda_envs_export
#
set -uo pipefail

INDIR="${1:-$HOME/conda_envs_export}"
LOGDIR="$INDIR/_logs"
FAILED="$INDIR/_FAILED.txt"
OK="$INDIR/_OK.txt"

if [ ! -d "$INDIR" ]; then
    echo "ERROR: no existe el directorio $INDIR"
    exit 1
fi

# Cargar conda en este shell
source "$HOME/miniforge3/etc/profile.d/conda.sh"

mkdir -p "$LOGDIR"
: > "$FAILED"
: > "$OK"

# mamba si está disponible, si no conda
SOLVER="mamba"
command -v mamba >/dev/null 2>&1 || SOLVER="conda"
echo ">>> Usando: $SOLVER"

shopt -s nullglob
for YML in "$INDIR"/*.yml; do
    ENV_NAME="$(basename "$YML" .yml)"

    if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
        echo "    [existe] $ENV_NAME  (saltando)"
        continue
    fi

    echo "    [crear] $ENV_NAME ..."
    if $SOLVER env create -f "$YML" -n "$ENV_NAME" > "$LOGDIR/${ENV_NAME}.log" 2>&1; then
        echo "        OK"
        echo "$ENV_NAME" >> "$OK"
    else
        echo "        FALLO  -> ver $LOGDIR/${ENV_NAME}.log"
        echo "$ENV_NAME" >> "$FAILED"
    fi
done

echo ""
echo "================================================================="
echo " Creados correctamente: $(wc -l < "$OK" | tr -d ' ')"
echo " Fallidos:              $(wc -l < "$FAILED" | tr -d ' ')"
[ -s "$FAILED" ] && { echo ""; echo " Revisar:"; cat "$FAILED"; }
echo ""
echo " Logs en: $LOGDIR"
echo "================================================================="