#!/usr/bin/env bash
#
# export_all_envs.sh  --  EJECUTAR EN TU MAC (local)
#
# Exporta todos tus entornos conda a ficheros YAML portables,
# aptos para recrearlos en Linux (farm22).
#
# Genera:
#   ./conda_envs_export/<nombre>.yml     (--from-history, portable)
#   ./conda_envs_export/_full/<nombre>.yml (export completo, de referencia)
#   ./conda_envs_export/_MANIFEST.txt
#
# Uso:
#   chmod +x export_all_envs.sh
#   ./export_all_envs.sh
#
set -uo pipefail

OUTDIR="./conda_envs_export"
FULLDIR="$OUTDIR/_full"
MANIFEST="$OUTDIR/_MANIFEST.txt"

# Entornos que NO se exportan (base y cualquiera que añadas aquí)
SKIP_ENVS=("base")

mkdir -p "$FULLDIR"
: > "$MANIFEST"

echo ">>> Detectando entornos..."

# Lista de rutas de entornos (columna final de 'conda env list')
# (compatible con el bash 3.2 que trae macOS: sin mapfile)
ENV_PATHS=()
while IFS= read -r line; do
    ENV_PATHS+=("$line")
done < <(conda env list | grep -v '^#' | grep -v '^[[:space:]]*$' | awk '{print $NF}')

for ENV_PATH in "${ENV_PATHS[@]}"; do
    ENV_NAME="$(basename "$ENV_PATH")"

    # Saltar entornos excluidos
    skip=0
    for s in "${SKIP_ENVS[@]}"; do
        [[ "$ENV_NAME" == "$s" ]] && skip=1
    done
    if [[ "$ENV_NAME" == "anaconda3" ]]; then skip=1; fi
    if [[ $skip -eq 1 ]]; then
        echo "    [skip] $ENV_NAME"
        continue
    fi

    echo "    [export] $ENV_NAME"

    # 1) Exportación PORTABLE (solo paquetes pedidos explícitamente).
    #    Es la que usarás para recrear en Linux.
    if conda env export --prefix "$ENV_PATH" --from-history > "$OUTDIR/${ENV_NAME}.yml" 2>/dev/null; then
        # Normalizar el campo 'name:' al nombre corto del entorno
        sed -i '' "1s|^name:.*|name: ${ENV_NAME}|" "$OUTDIR/${ENV_NAME}.yml" 2>/dev/null \
            || sed -i "1s|^name:.*|name: ${ENV_NAME}|" "$OUTDIR/${ENV_NAME}.yml"

        # Asegurar que bioconda/conda-forge estén como canales
        if ! grep -q "bioconda" "$OUTDIR/${ENV_NAME}.yml"; then
            python3 - "$OUTDIR/${ENV_NAME}.yml" <<'PY'
import sys, re
p = sys.argv[1]
txt = open(p).read()
if "channels:" in txt:
    txt = re.sub(r"channels:\n", "channels:\n  - conda-forge\n  - bioconda\n", txt, count=1)
else:
    txt = txt.replace("dependencies:", "channels:\n  - conda-forge\n  - bioconda\ndependencies:", 1)
open(p, "w").write(txt)
PY
        fi
        echo "$ENV_NAME" >> "$MANIFEST"
    else
        echo "        !! fallo exportando $ENV_NAME"
        continue
    fi

    # 2) Exportación COMPLETA (referencia: versiones exactas que usas en Mac)
    conda env export --prefix "$ENV_PATH" --no-builds > "$FULLDIR/${ENV_NAME}.yml" 2>/dev/null
done

echo ""
echo "================================================================="
echo " Exportados $(wc -l < "$MANIFEST" | tr -d ' ') entornos en: $OUTDIR"
echo ""
echo " Siguiente paso, copiar al farm:"
echo "   scp -r $OUTDIR jm79@farm22-head1:~/"
echo "================================================================="