#!/bin/bash
# ---------------------------------------------------------------------------
# srna_nivelB_collect.sh — FASE 2 / NIVEL B, paso 4: cierre, calibracion y cruce
#
#   ./srna_nivelB_collect.sh               # junta, clasifica, calibra, resume
#   ./srna_nivelB_collect.sh --objetivos lista.txt
#
# 1. Comprueba que estan todas las librerias (y dice que indices relanzar).
# 2. Junta las metricas por muestra en una tabla y la clasifica.
# 3. CALIBRACION con AEFE1: el control positivo de perfil EVE. Si AEFE1 no sale
#    EVE_NIRVS en la mayoria de las librerias, los umbrales de la clase EVE
#    estan mal y NO se interpreta nada mas (la leccion de dup_unmapped.sh:
#    primero una muestra de valor conocido).
# 4. Cruza con el RNA largo (small frente a largo, §4d del plan) usando la
#    cobertura por vOTU que ya existe.
# 5. Escribe RESUMEN.txt y las tablas que consumen las figuras.
#
# ANCLAJE: catalogo de 5,899 vOTUs (corte >=500 bp), decidido el 12 nov 2026.
# Las anotaciones que se usan son las que esa rama YA tiene; la tabla final de
# `_min500` todavia no existe (bloqueada por cov_fix), asi que es OPCIONAL: sin
# ella la clasificacion funciona igual, lo que falta son los nombres de virus.
#
# Variables (exportar antes):
#   QC_MIN500   por defecto $USER_ROOT/results/MERIDA/qc_control_min500
#   COV_LARGO   cobertura del RNA largo por vOTU (SIN encabezado)
#   VOTUS_RRNA  GFF de barrnap          ($R/rrna_min500/votus_rrna.gff)
#   VOTUS_GEN   cribado AaegL5          ($R/cribado_min500/votus_en_genoma.tsv)
#   TABLA       tabla final de la rama, si ya existe (opcional)
#   AEFE1_ID    nombre de la secuencia de AEFE1 en el indice (por defecto AEFE1)
#
# SALIDA  $RESULTS_DIR/nivelB/
#   nivelB_metricas.tsv      todas las muestras juntas
#   nivelB_clases.tsv        + clase y motivo por vOTU x muestra
#   nivelB_consenso.tsv      una fila por vOTU
#   nivelB_rpm_matriz.tsv    vOTU x muestra (para el bubble plot)
#   nivelB_small_vs_largo.tsv
#   RESUMEN.txt
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export S="${S:-$SCRIPT_DIR}"
source "$SCRIPT_DIR/srna_env.sh" "${SRNA_LOTE:-MERIDA_SRNA}" >/dev/null

OBJETIVOS=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --objetivos) OBJETIVOS="${2:?}"; shift 2 ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) echo "argumento desconocido: $1"; exit 2 ;;
  esac
done

R="$USER_ROOT/results/MERIDA"
QC_MIN500="${QC_MIN500:-$R/qc_control_min500}"
COV_LARGO="${COV_LARGO:-$QC_MIN500/cov/votu_coverage.tsv}"
VOTUS_RRNA="${VOTUS_RRNA:-$R/rrna_min500/votus_rrna.gff}"
VOTUS_GEN="${VOTUS_GEN:-$R/cribado_min500/votus_en_genoma.tsv}"
TABLA="${TABLA:-$QC_MIN500/TABLA_FINAL_votus_arb.tsv}"
AEFE1_ID="${AEFE1_ID:-AEFE1}"
LISTA="${SRNA_LIBS_FILE:-$SCRATCH/nivelB_libs.txt}"

OUT="$RESULTS_DIR/nivelB"
MET="$OUT/metricas"
[[ -d "$MET" ]] || { echo "ERROR: no existe $MET — corre antes srna_map_worker.sh"; exit 1; }
activate_env "$ENV_SRNA" >/dev/null 2>&1 || true

# --------------------------------------------------------------- 1. completitud
echo "== Nivel B: cierre de la fase 2"
faltan=(); ok=0
if [[ -s "$LISTA" ]]; then
  i=0
  while read -r lib; do
    [[ -n "$lib" ]] || continue
    i=$((i + 1))
    if [[ -s "$MET/$lib.votu_metrics.tsv" ]]; then ok=$((ok + 1)); else faltan+=("$i:$lib"); fi
  done < "$LISTA"
  echo "librerias con metricas: $ok/$i"
  if (( ${#faltan[@]} )); then
    echo "  !! faltan: ${faltan[*]}"
    echo "  relanzar solo esos indices:  bsub -J \"srnamap[$(IFS=,; echo "${faltan[*]%%:*}")]\" ..."
    exit 1
  fi
else
  echo "  !! sin lista $LISTA: se juntan los ficheros que haya en $MET"
fi

# --------------------------------------------------------------- 2. junta y clasifica
head -1 "$(ls "$MET"/*.votu_metrics.tsv | head -1)" > "$OUT/nivelB_metricas.tsv"
for f in "$MET"/*.votu_metrics.tsv; do tail -n +2 "$f"; done >> "$OUT/nivelB_metricas.tsv"
n_filas=$(( $(wc -l < "$OUT/nivelB_metricas.tsv") - 1 ))
echo "pares vOTU x muestra: $n_filas"

# listas de anotacion (opcionales, pero se avisa si faltan: cambian la lectura)
rrna_txt="$OUT/.votus_rrna.txt"; gen_txt="$OUT/.votus_en_genoma.txt"
if [[ -s "$VOTUS_RRNA" ]]; then
  awk -F'\t' '!/^#/ && NF > 3 {print $1}' "$VOTUS_RRNA" | sort -u > "$rrna_txt"
  echo "vOTUs con rRNA (barrnap): $(wc -l < "$rrna_txt")"
else
  : > "$rrna_txt"
  echo "  !! sin $VOTUS_RRNA: el 77-92 % de las bases mapeadas de esta rama va a"
  echo "     vOTUs de rRNA y SortMeRNA nunca corrio. Sin esta lista no se marcan."
fi
if [[ -s "$VOTUS_GEN" ]]; then
  awk 'NR > 1 {print $1}' "$VOTUS_GEN" | sort -u > "$gen_txt"
  echo "vOTUs en AaegL5 (cribado): $(wc -l < "$gen_txt")"
else
  : > "$gen_txt"; echo "  (sin $VOTUS_GEN: la columna en_cribado_genoma ira a 'no')"
fi

python3 "$SCRIPT_DIR/srna_clasificar.py" "$OUT/nivelB_metricas.tsv" \
  --salida "$OUT/nivelB_clases.tsv" --consenso "$OUT/nivelB_consenso.tsv" \
  --rrna "$rrna_txt" --en-genoma "$gen_txt"

# --------------------------------------------------------------- 3. calibracion AEFE1
cal=$(awk -F'\t' -v id="$AEFE1_ID" '
  NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "votu") cv = i; if ($i == "clase") cc = i }; next }
  $cv == id { n++; if ($cc == "EVE_NIRVS") e++ }
  END { printf "%d\t%d", n + 0, e + 0 }' "$OUT/nivelB_clases.tsv")
n_ae=${cal%%$'\t'*}; e_ae=${cal##*$'\t'}
if (( n_ae == 0 )); then
  cal_txt="SIN CALIBRAR — AEFE1 ('$AEFE1_ID') no tiene lecturas en ninguna libreria.
  Comprueba que esta en el indice viral (VIRAL_EXTRA) y que el ID coincide."
elif (( e_ae * 2 >= n_ae )); then
  cal_txt="PASA — AEFE1 sale EVE_NIRVS en $e_ae de $n_ae librerias."
else
  cal_txt="NO PASA — AEFE1 sale EVE_NIRVS en solo $e_ae de $n_ae librerias.
  Los umbrales de la clase EVE estan mal calibrados: NO interpretar el resto
  hasta revisarlos (mira sus filas en nivelB_clases.tsv y el motivo)."
fi
echo; echo "calibracion AEFE1: ${cal_txt%%$'\n'*}"

# --------------------------------------------------------------- 4. matriz y small vs largo
python3 - "$OUT/nivelB_metricas.tsv" "$OUT/nivelB_rpm_matriz.tsv" <<'PY'
import csv, sys
from collections import defaultdict
rpm = defaultdict(dict); muestras = set()
with open(sys.argv[1]) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        rpm[r["votu"]][r["muestra"]] = r["rpm"]; muestras.add(r["muestra"])
ms = sorted(muestras)
with open(sys.argv[2], "w") as out:
    out.write("votu\t" + "\t".join(ms) + "\n")
    for v in sorted(rpm):
        out.write(v + "\t" + "\t".join(rpm[v].get(m, "0") for m in ms) + "\n")
PY
echo "matriz RPM: $(( $(wc -l < "$OUT/nivelB_rpm_matriz.tsv") - 1 )) vOTUs"

# votu_coverage.tsv NO tiene encabezado (trampa conocida):
#   sample, votu, longitud, bases_cubiertas, pct_cobertura, profundidad_media
# El emparejamiento libreria de sRNA -> muestra de RNA largo sale del samplesheet
# (merida_id), que es la unica fuente de verdad: el indice S## no coincide.
if [[ -s "$COV_LARGO" && -s "$SRNA_SHEET" ]]; then
  awk -F'\t' '
    NR == FNR { if (FNR > 1 && $2 != "" && $2 != "NA") par[$1] = $2; next }
    FNR == 1 && FILENAME ~ /nivelB_metricas/ { for (i = 1; i <= NF; i++) h[$i] = i; next }
    FILENAME ~ /nivelB_metricas/ { m = $(h["muestra"]); if (m in par) rpm[par[m] "\t" $(h["votu"])] = $(h["rpm"]); next }
    { k = $1 "\t" $2; if (k in rpm) print k "\t" rpm[k] "\t" $5 "\t" $6 }
  ' "$SRNA_SHEET" "$OUT/nivelB_metricas.tsv" "$COV_LARGO" \
  | sort -k1,1 -k2,2 \
  | cat <(printf 'muestra_largo\tvotu\trpm_srna\tpct_cobertura_largo\tprofundidad_largo\n') - \
  > "$OUT/nivelB_small_vs_largo.tsv"
  echo "small vs largo: $(( $(wc -l < "$OUT/nivelB_small_vs_largo.tsv") - 1 )) pares emparejados"
else
  echo "  (sin $COV_LARGO o sin samplesheet: no se hace el cruce small vs largo)"
fi

# --------------------------------------------------------------- 5. nombres de virus si los hay
if [[ -s "$TABLA" ]]; then
  awk -F'\t' '
    NR == FNR { if (FNR == 1) { for (i = 1; i <= NF; i++) { if ($i == "votu") cv = i
                                                            if ($i == "nombre_virus") cn = i
                                                            if ($i == "especie" && !cn) cn = i } ; next }
                if (cn && $cn != "" && $cn != "NA") nom[$cv] = $cn; next }
    FNR == 1 { print $0 "\tnombre_virus"; next }
    { print $0 "\t" ($1 in nom ? nom[$1] : "NA") }
  ' "$TABLA" "$OUT/nivelB_consenso.tsv" > "$OUT/.consenso.tmp" \
    && mv "$OUT/.consenso.tmp" "$OUT/nivelB_consenso.tsv"
  echo "consenso anotado con nombre_virus desde $(basename "$TABLA")"
else
  echo "  (sin tabla final de la rama _min500: el consenso va sin nombre_virus)"
fi

# --------------------------------------------------------------- 6. RESUMEN
{
  echo "NIVEL B — cruce por secuencia del ARN pequeno contra las vOTUs"
  echo "fecha: $(date '+%F %T')"
  echo "catalogo: $(basename "$(dirname "$BT1_VIRAL")") (anclaje: 5,899 vOTUs, corte >=500 bp)"
  echo "librerias: $ok   pares vOTU x muestra: $n_filas"
  echo
  echo "CALIBRACION"
  echo "  $cal_txt"
  echo
  echo "DENOMINADOR DE RPM por libreria (18-35 nt no rRNA)"
  awk -F'\t' '$2 == "denominador_rpm" {printf "  %-6s %12d\n", $1, $3}' "$OUT"/particion/*.particion.tsv | sort
  echo
  echo "rRNA fuera, por libreria"
  awk -F'\t' '$2 == "lecturas_18_35" {t[$1] = $3} $2 == "rrna" {r[$1] = $3}
       END {for (l in t) printf "  %-6s %5.1f %%\n", l, 100 * r[l] / t[l]}' \
    "$OUT"/particion/*.particion.tsv | sort
  echo
  echo "CONSENSO por vOTU"
  awk -F'\t' 'NR > 1 {n[$10]++} END {for (c in n) printf "  %-24s %5d\n", c, n[c]}' \
    "$OUT/nivelB_consenso.tsv" | sort -k2,2nr
  echo
  echo "Los MIXTO_replica_y_EVE son los que importan: replica en unos pools y perfil"
  echo "EVE en otros. Es la firma que se esperaba de CFAV (infeccion real en M72,"
  echo "M61 y M41 sobre un fondo de ~27 % de cobertura en todas las librerias)."
  echo
  echo "vOTUs con mas RPM (top 30, sin rRNA)"
  awk -F'\t' 'NR > 1 && $10 != "no_interpretable_rrna" {printf "  %-18s %10.1f  %-24s %s\n", $1, $9, $10, (NF >= 11 ? $11 : "")}' \
    "$OUT/nivelB_consenso.tsv" | sort -k2,2nr | head -30
  if [[ -n "$OBJETIVOS" && -s "$OBJETIVOS" ]]; then
    echo
    echo "OBJETIVOS ($(basename "$OBJETIVOS")) — las preguntas que el nivel A no pudo cerrar"
    awk -F'\t' 'NR == FNR {q[$1] = 1; next}
      FNR == 1 {next}
      ($1 in q) {printf "  %-18s muestras %-3s replica %-3s EVE %-3s amb %-3s rpm_max %9.1f  %s\n", $1, $2, $3, $5, $6, $9, $10}' \
      "$OBJETIVOS" "$OUT/nivelB_consenso.tsv"
    echo "  sin lecturas en ninguna libreria:"
    awk -F'\t' 'NR == FNR { if (FNR > 1) vista[$1] = 1; next }
                $1 != "" && $1 != "votu" && !($1 in vista) {print "    " $1}' \
      "$OUT/nivelB_consenso.tsv" "$OBJETIVOS"
  fi
  echo
  echo "AVISOS QUE SIGUEN VIGENTES"
  echo "  1. El 45 % de las vOTUs de esta rama cae en AaegL5 y el catalogo es ~87 % huesped:"
  echo "     solo se interpretan el subconjunto viral y los candidatos a EVE/NIRVS."
  echo "  2. El RPM va sobre 18-35 no rRNA y depende de cuanto miRNA tenga cada libreria."
  echo "     Reportar tambien normalizado por miRNA; si las dos ordenan distinto, decirlo."
  echo "  3. frac_unica baja = la lectura la recogen varias vOTUs emparentadas."
  echo "     Contar taxones sin mirar esa columna infla el lado del ARN pequeno."
  echo "  4. Un EVE en region repetida tiene lecturas en decenas de loci:"
  echo "     loci_genoma_mediana y frac_en_genoma lo hacen visible, no se esconde."
} > "$OUT/RESUMEN.txt"

rm -f "$rrna_txt" "$gen_txt"
echo
sed -n '1,/^CONSENSO/p' "$OUT/RESUMEN.txt"
echo "..."
echo "Done: $OUT/RESUMEN.txt ($((SECONDS/60)) min)"
