#!/bin/bash
# ---------------------------------------------------------------------------
# srna_b0b_virus.sh — intercambio M41 <-> M48 con VARIANTES VIRALES (B0b)
#
# B0 (ADNmt) no se pudo calibrar: la variacion mitocondrial entre pools es baja y hay sesgos
# propios de cada tecnica (modificaciones de ARNt, piRNA de NUMTs). Los virus son mucho mas
# especificos de cada pool: cada pool tiene su propia version de Verdadero, Chaq-like, etc.
#
#  1. En la union de contigs del RNA total se buscan, para cada pool, sus contigs de cada virus
#     marcador (trozos de 150 nt alineados contra la vOTU del virus).
#  2. Las lecturas de ARN pequeno de cada libreria se alinean contra todos esos contigs a la vez.
#  3. Una lectura cuenta solo si su region esta ensamblada en TODOS los pools y coincide mejor
#     con UNO solo. Se asigna a ese pool.
#  4. Calibracion: las librerias de control deben asignar la mayoria a su propio pool.
#
# Uso (desde la carpeta de scripts, sin source):
#   ./srna_b0b_virus.sh check
#   bsub -J b0b -o <logs>/b0b.%J.log -e <logs>/b0b.%J.err -q normal -n 4 -M 8000 \
#        -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" "$PWD/srna_b0b_virus.sh run"
#
# Variables (exportar antes para cambiar):
#   SRNA_ROOT  UNION_FA  TABLA  VOTUS_FA  SRNA_LIBS  RNA_POOLS  PAREJAS
#   MARCADORES  virus separados por ';' tal como aparecen en la columna nombre_virus de TABLA
#
# SALIDA: $RESULTS_DIR/b0b_virus/  (RESULTS_DIR de MERIDA_SRNA)
#   RESUMEN.txt  calibracion.tsv  fracciones.tsv  asignaciones.tsv  marcadores_por_pool.tsv
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export S="${S:-$SCRIPT_DIR}"
source "$SCRIPT_DIR/srna_env.sh" MERIDA_SRNA >/dev/null

cmd="${1:?uso: srna_b0b_virus.sh check|run}"
THREADS="${LSB_DJOB_NUMPROC:-2}"
R="$USER_ROOT/results/MERIDA"

SRNA_ROOT="${SRNA_ROOT:-/lustre/scratch126/tol/teams/lawniczak/projects/vector_eve/data/MosquitoSmallRNA}"
UNION_FA="${UNION_FA:-$R/contig_annot/contigs_union_MERIDA.fasta}"
VOTUS_FA="${VOTUS_FA:-$REFS_DIR/votu_index/votus_from_index.fasta}"
SRNA_LIBS="${SRNA_LIBS:-M41 M48 M2 M11 M21 M60 M72}"
RNA_POOLS="${RNA_POOLS:-M48_S5 M2_S2 M11_S3 M21_S29 M60_S11 M72_S20}"
PAREJAS="${PAREJAS:-M2=M2_S2,M11=M11_S3,M21=M21_S29,M60=M60_S11,M72=M72_S20,M41=M48_S5}"
MARCADORES="${MARCADORES:-Verdadero virus;Chaq-like virus;Guadeloupe mosquito virus;Aedes anphevirus}"
OUT="$RESULTS_DIR/b0b_virus"
WORK="$SCRATCH/b0b_virus"

# TABLA: la tabla final arbitrada si existe, si no la tabla final
if [[ -z "${TABLA:-}" ]]; then
  for t in "$R/qc_control/TABLA_FINAL_votus_arb.tsv" "$R/TABLA_FINAL_votus_arb.tsv" \
           "$R/qc_control/TABLA_FINAL_votus.tsv" "$R/TABLA_FINAL_votus.tsv"; do
    [[ -s "$t" ]] && { TABLA="$t"; break; }
  done
fi
TABLA="${TABLA:-}"

# --------------------------------------------------------------- utilidades (iguales que en B0)
find_srna() { ls "$SRNA_ROOT"/*/00_CleanData/"$1".clean.fa.gz 2>/dev/null | head -1 || true; }
formato_fa() {
  local h; h=$( { zcat "$1" 2>/dev/null || true; } | head -1)
  if   [[ "$h" =~ _x[0-9]+$ ]];          then echo "colapsado_x"
  elif [[ "$h" =~ [[:space:]][0-9]+$ ]]; then echo "colapsado_espacio"
  else echo "expandido"; fi
}
lecturas_18_35() {
  local fa="$1" modo; modo=$(formato_fa "$fa")
  zcat "$fa" | awk -v modo="$modo" '
    /^>/ { n = 1
           if (modo == "colapsado_x")            { if (match($0, /_x[0-9]+$/)) n = substr($0, RSTART + 2) + 0 }
           else if (modo == "colapsado_espacio") { n = $NF + 0 }
           next }
    { L = length($0); if (L >= 18 && L <= 35) for (i = 0; i < n; i++) { k++; print ">r" k; print $0 } }'
}

# vOTUs de cada virus marcador segun TABLA: imprime "votu<TAB>virus_sin_espacios"
votus_marcadores() {
  awk -F'\t' -v lista="$MARCADORES" '
    NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "votu") cv = i; if ($i == "nombre_virus") cn = i; if ($i == "especie") ce = i }
              if (!cn) cn = ce
              n = split(lista, L, ";"); for (k = 1; k <= n; k++) want[L[k]] = 1; next }
    ($cn in want) { v = $cn; gsub(/ /, "_", v); print $cv "\t" v }' "$TABLA"
}

# --------------------------------------------------------------- check
if [[ "$cmd" == check ]]; then
  echo "== B0b: intercambio M41/M48 con variantes virales — comprobacion previa"
  echo
  activate_env "$ENV_SRNA" >/dev/null 2>&1 || echo "  !! no se activa el entorno '$ENV_SRNA'"
  for t in bowtie bowtie-build samtools python3; do
    printf '  %-13s %s\n' "$t" "$(command -v "$t" >/dev/null && "$t" --version 2>&1 | head -1 || echo 'NO ESTA')"
  done
  printf '  %-13s %s\n' seqkit "$(command -v seqkit >/dev/null && seqkit version 2>&1 | head -1 || echo 'NO ESTA')"
  for f in srna_b0b_asignar.py srna_b0b_resumen.py; do
    [[ -s "$SCRIPT_DIR/$f" ]] && echo "  $f ok" || echo "  !! falta $SCRIPT_DIR/$f"
  done
  echo
  echo "-- union de contigs: $UNION_FA"
  if [[ -s "$UNION_FA" ]]; then
    for p in $RNA_POOLS; do printf '  %-8s %6d contigs\n' "$p" "$(grep -c "^>${p}__" "$UNION_FA" || true)"; done
  else echo "  !! no existe (exporta UNION_FA=<ruta>)"; fi
  echo
  echo "-- tabla de vOTUs: ${TABLA:-NO ENCONTRADA (exporta TABLA=<ruta>)}"
  echo "-- secuencias de vOTUs: $VOTUS_FA $([[ -s "$VOTUS_FA" ]] && echo ok || echo '!! NO EXISTE')"
  if [[ -n "$TABLA" && -s "$TABLA" ]]; then
    echo "   vOTUs por virus marcador:"
    votus_marcadores | cut -f2 | sort | uniq -c | sed 's/^/     /'
    falta=$(IFS=';'; for v in $MARCADORES; do vv=${v// /_}; votus_marcadores | cut -f2 | grep -qx "$vv" || echo "$v"; done)
    [[ -z "$falta" ]] || echo "   !! sin vOTU en la tabla: $falta"
  fi
  echo
  echo "-- ARN pequeno"
  for l in $SRNA_LIBS; do
    f=$(find_srna "$l")
    [[ -n "$f" ]] && printf '  %-5s %s (%s)\n' "$l" "${f#$SRNA_ROOT/}" "$(formato_fa "$f")" || printf '  %-5s !! NO ENCONTRADO\n' "$l"
  done
  exit 0
fi
[[ "$cmd" == run ]] || { echo "subcomando desconocido: $cmd"; exit 2; }

mkdir -p "$OUT" "$WORK"
echo "B0b virus — $(date '+%F %T') en $(hostname), $THREADS hilos"
[[ -s "$UNION_FA" && -n "$TABLA" && -s "$TABLA" && -s "$VOTUS_FA" ]] || { echo "ERROR: faltan entradas; corre 'check'"; exit 1; }
activate_env "$ENV_SRNA"

# --------------------------------------------------------------- 1. referencias de los virus marcadores
votus_marcadores > "$WORK/votus_marcadores.tsv"
[[ -s "$WORK/votus_marcadores.tsv" ]] || { echo "ERROR: ningun virus marcador en $TABLA"; exit 1; }
seqkit grep -f <(cut -f1 "$WORK/votus_marcadores.tsv") "$VOTUS_FA" \
  | awk 'NR==FNR{v[$1]=$2; next} /^>/{id=substr($1,2); print ">" v[id] "__" id; next} {print}' "$WORK/votus_marcadores.tsv" - \
  > "$WORK/virus_refs.fa"
echo "referencias: $(grep -c '^>' "$WORK/virus_refs.fa") vOTUs de $(cut -f2 "$WORK/votus_marcadores.tsv" | sort -u | wc -l) virus"
bowtie-build -q "$WORK/virus_refs.fa" "$WORK/virus_refs"

# --------------------------------------------------------------- 2. contigs de cada virus en cada pool
patron=$(echo $RNA_POOLS | sed 's/ /|/g')
seqkit grep -r -p "^(${patron})__" "$UNION_FA" > "$WORK/union_pools.fa"
# trozos de 150 nt cada 75; un contig es del virus si >= 50 % de sus trozos alinean (<= 3 errores)
seqkit sliding -W 150 -s 75 "$WORK/union_pools.fa" 2>/dev/null | seqkit seq -m 150 > "$WORK/trozos.fa"
bowtie -p "$THREADS" -f -v 3 -k 1 --best -S --no-unal "$WORK/virus_refs" "$WORK/trozos.fa" 2> "$WORK/trozos.log" \
  | awk -F'\t' '!/^@/ {c = $1; sub(/_sliding:.*/, "", c); split($3, a, "__"); n[c "\t" a[1]]++} END{for (k in n) print k "\t" n[k]}' \
  > "$WORK/trozos_por_contig.tsv"
grep '^>' "$WORK/trozos.fa" | sed 's/^>//; s/_sliding:.*//; s/[[:space:]].*//' | sort | uniq -c | awk '{print $2 "\t" $1}' > "$WORK/trozos_totales.tsv"
awk -F'\t' 'NR==FNR{tot[$1]=$2; next} {if ($3 / tot[$1] >= 0.5) print $1 "\t" $2}' \
  "$WORK/trozos_totales.tsv" "$WORK/trozos_por_contig.tsv" \
  | sort -k1,1 -k3,3nr | awk '!vista[$1]++' > "$WORK/contig_virus.tsv"
# FASTA de marcadores con cabecera pool__virus__contig
awk -F'\t' 'NR==FNR{v[$1]=$2; next}
  /^>/ {id = substr($1, 2); keep = (id in v); if (keep) {split(id, a, "__"); print ">" a[1] "__" v[id] "__" id}; next}
  keep' "$WORK/contig_virus.tsv" "$WORK/union_pools.fa" > "$WORK/marcadores.fa"
awk -F'\t' '{split($1, a, "__"); print $2 "\t" a[1]}' "$WORK/contig_virus.tsv" | sort -u > "$WORK/pools_por_virus.tsv"
{ printf 'virus\tpool\tcontigs\tlongitud_total\n'
  seqkit fx2tab -n -l "$WORK/marcadores.fa" | awk -F'\t' '{split($1, a, "__"); k = a[2] "\t" a[1]; n[k]++; L[k] += $NF} END{for (k in n) print k "\t" n[k] "\t" L[k]}' | sort
} > "$OUT/marcadores_por_pool.tsv"
echo; echo "contigs marcadores por virus y pool:"; awk -F"\t" '{printf "  %-28s %-9s %8s %14s\n", $1, $2, $3, $4}' "$OUT/marcadores_por_pool.tsv"
bowtie-build -q "$WORK/marcadores.fa" "$WORK/marcadores"

# --------------------------------------------------------------- 3. lecturas de ARN pequeno
printf 'libreria\tvirus\tpool\tlecturas\n' > "$OUT/asignaciones.tsv"
echo
for l in $SRNA_LIBS; do
  f=$(find_srna "$l"); [[ -n "$f" ]] || { echo "ERROR: no encuentro $l"; exit 1; }
  lecturas_18_35 "$f" > "$WORK/$l.18_35.fa"
  # -v 2 -a: todas las alineaciones con <= 2 errores. Sin --reorder: con bowtie 1.3.1 y -p se cuelga.
  # bowtie escribe juntas las alineaciones de cada lectura; el script de asignacion lo verifica.
  bowtie -p "$THREADS" -f -v 2 -a -S --no-unal "$WORK/marcadores" "$WORK/$l.18_35.fa" 2> "$WORK/$l.bowtie.log" \
    | python3 "$SCRIPT_DIR/srna_b0b_asignar.py" "$l" "$WORK/pools_por_virus.tsv" >> "$OUT/asignaciones.tsv"
  rm -f "$WORK/$l.18_35.fa"
done

# --------------------------------------------------------------- 4. resumen y calibracion
echo
python3 "$SCRIPT_DIR/srna_b0b_resumen.py" "$OUT/asignaciones.tsv" "$OUT" --parejas "$PAREJAS" | tee "$OUT/RESUMEN.txt"
cat >> "$OUT/RESUMEN.txt" <<EOF

Lectura: cada libreria de control debe asignar la mayor fraccion de sus lecturas informativas a su
propio pool. Si es asi y "M41" asigna la mayor fraccion a M48_S5, el intercambio queda confirmado
por secuencia. "M48" no tiene virus (0 contigs virales en el reporte): se espera n muy bajo.
EOF
echo
echo "Done: $OUT ($((SECONDS/60)) min)"
