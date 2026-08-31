#!/bin/bash
# ---------------------------------------------------------------------------
# clstr_aporte.sh — que aporta REALMENTE cada ensamblador, leyendo los .clstr
#                   que deja derep_contigs.sh.
#
# POR QUE NO VALE CONTAR REPRESENTANTES
#   CD-HIT elige como representante la secuencia mas larga del grupo, y en los
#   EMPATES gana la que entro primero en el fichero. derep_contigs.sh concatena
#   siempre en el mismo orden, con rnaviral el primero, asi que su recuento de
#   representantes esta inflado por construccion. No mide aportacion.
#
# LO QUE SI MIDE
#   SOLO_EL   grupos donde ese ensamblador es el UNICO presente. Eso es material
#             que nadie mas encontro.
#   GANA      grupos donde es el representante Y hay mas de un ensamblador, o sea
#             donde dio la version mas larga de algo que otros tambien vieron.
#             Esta es la columna que revela quien ensambla mas completo.
#
# Uso:
#   ./clstr_aporte.sh                 # todas las muestras
#   ./clstr_aporte.sh 3000            # solo grupos con representante >=3000 bp
#   ./clstr_aporte.sh 3000 M71_S19    # una muestra
# ---------------------------------------------------------------------------
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config.sh"

MINLEN="${1:-0}"
MUESTRA="${2:-}"
UNION="${UNION_DIR:-$SCRATCH/union}"

if [[ -n "$MUESTRA" ]]; then
  files=( "$UNION/$MUESTRA/$MUESTRA.clstr" )
else
  mapfile -t files < <(find "$UNION" -name "*.clstr" -size +0 | sort)
fi
(( ${#files[@]} )) || { echo "No hay ficheros .clstr en $UNION"; exit 1; }

echo "==========================================================="
echo "clstr_aporte.sh   ficheros: ${#files[@]}   longitud minima: $MINLEN bp"
echo "==========================================================="

awk -v MIN="$MINLEN" '
function flush(   k, c, only) {
  if (!nmiem) return
  if (replen < MIN) { return }
  grupos++
  c = 0
  for (k in tags) { c++; only = k }
  for (k in tags) presente[k]++
  if (c == 1) { solo[only]++; if (replen > solomax[only]) solomax[only] = replen }
  else if (rep != "") { gana[rep]++ }
}
/^>Cluster/ { flush(); delete tags; nmiem = 0; rep = ""; replen = 0; next }
{
  if (!match($0, />[A-Za-z0-9-]+__/)) next
  t = substr($0, RSTART + 1, RLENGTH - 3)
  tags[t] = 1
  nmiem++
  L = 0
  if (match($0, /[0-9]+nt,/)) L = substr($0, RSTART, RLENGTH - 3) + 0
  if ($0 ~ /\*[[:space:]]*$/) { rep = t; replen = L }
}
END {
  flush()
  printf "%-12s %10s %10s %10s %12s\n", "ensamblador", "en_grupos", "SOLO_EL", "GANA", "solo_max_bp"
  for (k in presente)
    printf "%-12s %10d %10d %10d %12d\n", k, presente[k], solo[k]+0, gana[k]+0, solomax[k]+0
  printf "\ngrupos contados: %d\n", grupos
}' "${files[@]}" | { read -r h; echo "  $h"; sort -k3 -rn | sed 's/^/  /'; }

cat <<'EOF'

  SOLO_EL      grupos donde ese ensamblador es el unico presente
               -> material que nadie mas encontro
  GANA         grupos compartidos donde ese ensamblador dio la version mas larga
               -> quien ensambla mas completo
  solo_max_bp  el contig mas largo entre los que solo el encontro

  Los representantes a secas NO se cuentan aqui: en los empates gana el primero
  del fichero, que siempre es rnaviral, y eso inflaria su cifra.
EOF
