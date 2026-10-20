#!/bin/bash
# ---------------------------------------------------------------------------
# srna_b0_intercambio.sh — confirma o descarta el intercambio de etiquetas M41 <-> M48
# en las librerias de ARN pequeno, con frecuencias alelicas del genoma mitocondrial.
#
# Idea: cada pool de 15 mosquitos tiene su propia mezcla de haplotipos mitocondriales.
# Las lecturas de ARN pequeno y las de RNA total del MISMO pool deben dar las mismas
# frecuencias en los sitios variables del ADNmt. Se calibra con 5 parejas que no estan
# en duda (M2, M11, M21, M60, M72): cada una debe parecerse sobre todo a su propio pool.
# Si "M41" se parece sobre todo a M48_S5, el intercambio queda confirmado.
#
# Uso:
#   ./srna_b0_intercambio.sh check      # nodo de login, ~1 min: herramientas, entradas, formatos
#   bsub -J b0 -o "$LOGS_DIR/b0.%J.log" -e "$LOGS_DIR/b0.%J.err" \
#        -q normal -n 4 -M 8000 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "$S/srna_b0_intercambio.sh run"
#
# Variables que se pueden exportar antes:
#   SRNA_ROOT       carpeta con Mexico1/ Mexico2/ ... (00_CleanData/<lib>.clean.fa.gz)
#   SRNA_LIBS       librerias de ARN pequeno     (por defecto: M41 M48 M2 M11 M21 M60 M72)
#   RNA_POOLS       pools de RNA total           (por defecto: M48_S5 M2_S2 M11_S3 M21_S29 M60_S11 M72_S20)
#   RNASEQ_BAM_DIR  carpeta de los BAM de bowtie2 contra AaegL5 (si no, se busca en las rutas habituales)
#   MITO_ID         nombre del cromosoma mitocondrial en el indice (si no, se detecta)
#
# SALIDA: $RESULTS_DIR/b0_intercambio/  (RESULTS_DIR de MERIDA_SRNA)
#   emparejamiento.tsv  distancias.tsv  frecuencias_sitios_variables.tsv  conteos.tsv  RESUMEN.txt
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export S="${S:-$SCRIPT_DIR}"
source "$SCRIPT_DIR/srna_env.sh" MERIDA_SRNA >/dev/null

cmd="${1:?uso: srna_b0_intercambio.sh check|run}"
THREADS="${LSB_DJOB_NUMPROC:-2}"

SRNA_ROOT="${SRNA_ROOT:-/lustre/scratch126/tol/teams/lawniczak/projects/vector_eve/data/MosquitoSmallRNA}"
SRNA_LIBS="${SRNA_LIBS:-M41 M48 M2 M11 M21 M60 M72}"
RNA_POOLS="${RNA_POOLS:-M48_S5 M2_S2 M11_S3 M21_S29 M60_S11 M72_S20}"
# parejas esperadas: controles primero, la hipotesis al final
PAREJAS="${PAREJAS:-M2=M2_S2,M11=M11_S3,M21=M21_S29,M60=M60_S11,M72=M72_S20,M41=M48_S5}"
OUT="$RESULTS_DIR/b0_intercambio"
WORK="$SCRATCH/b0_intercambio"
MITO_FA="$SRNA_REFS/fasta/mito_AaegL5.fa"
BT1_MITO="$SRNA_REFS/bt1_mito/mito"
PY="$SCRIPT_DIR/srna_b0_frecuencias.py"

# --------------------------------------------------------------- utilidades
find_srna() { ls "$SRNA_ROOT"/*/00_CleanData/"$1".clean.fa.gz 2>/dev/null | head -1 || true; }

find_bam() {
  local d f
  for d in ${RNASEQ_BAM_DIR:-} "$USER_ROOT/work/MERIDA/mapped/MERIDA" "$USER_ROOT/results/MERIDA/mapped_bam/MERIDA" \
           "$USER_ROOT/work/MERIDA/mapped_bam/MERIDA" "$USER_ROOT/work/MERIDA/mapped"; do
    [[ -d "$d" ]] || continue
    for f in "$d/$1.bam" "$d/$1".*.bam "$d/$1"_*.bam; do
      [[ -s "$f" ]] && { echo "$f"; return 0; }
    done
  done
  return 0
}

# Formato del FASTA limpio: "expandido" (una secuencia por lectura) o "colapsado"
# (secuencias unicas con su conteo en la cabecera: >xxx_x123  o  >xxx 123).
formato_fa() {
  local h; h=$( { zcat "$1" 2>/dev/null || true; } | head -1)
  if   [[ "$h" =~ _x[0-9]+$ ]];              then echo "colapsado_x"
  elif [[ "$h" =~ [[:space:]][0-9]+$ ]];     then echo "colapsado_espacio"
  else echo "expandido"; fi
}

# Lecturas 18-35 nt, expandidas segun su conteo, a stdout en FASTA.
lecturas_18_35() {
  local fa="$1" modo; modo=$(formato_fa "$fa")
  zcat "$fa" | awk -v modo="$modo" '
    /^>/ { n = 1
           if (modo == "colapsado_x")       { if (match($0, /_x[0-9]+$/)) n = substr($0, RSTART + 2) + 0 }
           else if (modo == "colapsado_espacio") { n = $NF + 0 }
           next }
    { L = length($0); if (L >= 18 && L <= 35) for (i = 0; i < n; i++) { k++; print ">r" k; print $0 } }'
}

detect_mito_id() {
  [[ -n "${MITO_ID:-}" ]] && { echo "$MITO_ID"; return; }
  { bowtie2-inspect -n "$BT2_INDEX" 2>/dev/null || true; } | grep -iE 'NC_035159|mitochond|^MT([[:space:]]|$)|^chrM' | awk 'NR==1{print $1}' || true
}

bam_indexado() { [[ -s "$1.bai" || -s "${1%.bam}.bai" || -s "$1.csi" ]]; }

# --------------------------------------------------------------- check
if [[ "$cmd" == check ]]; then
  echo "== B0: intercambio M41/M48 — comprobacion previa"
  echo
  echo "-- herramientas"
  activate_env "$ENV_SRNA" >/dev/null 2>&1 || echo "  !! no se activa el entorno '$ENV_SRNA' (crealo antes, ver PLAN_smallRNA.md)"
  for t in bowtie bowtie-build samtools python3; do
    printf '  %-13s %s\n' "$t" "$(command -v "$t" >/dev/null && "$t" --version 2>&1 | head -1 || echo 'NO ESTA')"
  done
  [[ -s "$PY" ]] && echo "  srna_b0_frecuencias.py  ok" || echo "  !! falta $PY (copialo junto a este script)"
  echo
  echo "-- ARN pequeno ($SRNA_ROOT)"
  for l in $SRNA_LIBS; do
    f=$(find_srna "$l")
    if [[ -n "$f" ]]; then
      printf '  %-5s %-60s %8s  formato: %s\n' "$l" "${f#$SRNA_ROOT/}" "$(du -h "$f" | cut -f1)" "$(formato_fa "$f")"
    else
      printf '  %-5s !! NO ENCONTRADO\n' "$l"
    fi
  done
  l1=$(echo $SRNA_LIBS | awk '{print $1}'); f=$(find_srna "$l1")
  [[ -n "$f" ]] && { echo "  primeras lineas de $l1:"; { zcat "$f" 2>/dev/null || true; } | head -4 | sed 's/^/      /'; }
  echo
  echo "-- RNA total (BAM de bowtie2 contra AaegL5)"
  activate_env "$ENV_BOWTIE2" >/dev/null 2>&1 || echo "  !! no se activa $ENV_BOWTIE2"
  for p in $RNA_POOLS; do
    b=$(find_bam "$p")
    if [[ -n "$b" ]]; then
      so=$(samtools view -H "$b" 2>/dev/null | awk '/^@HD/{for(i=1;i<=NF;i++) if($i~/^SO:/) print substr($i,4)}')
      printf '  %-8s %-70s %8s  orden: %-10s indice: %s\n' "$p" "$b" "$(du -h "$b" | cut -f1)" "${so:-?}" \
             "$(bam_indexado "$b" && echo si || echo no)"
    else
      printf '  %-8s !! NO ENCONTRADO (exporta RNASEQ_BAM_DIR=<carpeta>)\n' "$p"
    fi
  done
  echo
  echo "-- cromosoma mitocondrial en $BT2_INDEX"
  mid=$(detect_mito_id)
  if [[ -n "$mid" ]]; then echo "  MITO_ID = $mid"
  else echo "  !! no lo detecto; nombres que parecen candidatos:"
       { bowtie2-inspect -n "$BT2_INDEX" 2>/dev/null || true; } | awk 'length($1)<20' | head -10 | sed 's/^/      /'
       echo "     exporta MITO_ID=<nombre> y repite"; fi
  b=$(find_bam "$(echo $RNA_POOLS | awk '{print $1}')")
  if [[ -n "$b" && -n "$mid" ]]; then
    hdr=$(samtools view -H "$b" 2>/dev/null || true)
    grep -q "SN:$mid[[:space:]]" <<< "$hdr" && echo "  el BAM de RNA total tiene la secuencia $mid: ok" \
      || echo "  !! el BAM de RNA total no tiene una secuencia llamada $mid"
  fi
  [[ -s "$MITO_FA" ]] && echo "  FASTA mitocondrial ya extraido: $MITO_FA" || echo "  FASTA mitocondrial: se extraera en 'run'"
  exit 0
fi

[[ "$cmd" == run ]] || { echo "subcomando desconocido: $cmd"; exit 2; }

mkdir -p "$OUT" "$WORK"
echo "B0 intercambio — $(date '+%F %T') en $(hostname), $THREADS hilos"

# --------------------------------------------------------------- 1. referencia mitocondrial
activate_env "$ENV_BOWTIE2"
MITO_ID=$(detect_mito_id)
[[ -n "$MITO_ID" ]] || { echo "ERROR: no detecto el cromosoma mitocondrial; exporta MITO_ID"; exit 1; }
echo "cromosoma mitocondrial: $MITO_ID"
if [[ ! -s "$MITO_FA" ]]; then
  mkdir -p "$(dirname "$MITO_FA")"
  echo "extrayendo $MITO_ID de $BT2_INDEX (recorre el genoma entero, unos minutos)"
  bowtie2-inspect "$BT2_INDEX" | awk -v id="$MITO_ID" '/^>/{p = ($1 == ">" id)} p' > "$MITO_FA.tmp"
  mv "$MITO_FA.tmp" "$MITO_FA"
fi
L=$(awk '!/^>/{n += length($0)} END{print n}' "$MITO_FA")
echo "  $MITO_FA: $L pb"
(( L > 14000 && L < 20000 )) || { echo "ERROR: longitud mitocondrial rara ($L pb)"; exit 1; }
samtools faidx "$MITO_FA"

: > "$WORK/nombres.txt"
printf 'muestra\ttipo\tfuente\tlecturas_entrada\tlecturas_mito\n' > "$OUT/conteos.tsv"

# --------------------------------------------------------------- 2. RNA total: lecturas del ADNmt
for p in $RNA_POOLS; do
  b=$(find_bam "$p"); [[ -n "$b" ]] || { echo "ERROR: no encuentro el BAM de $p"; exit 1; }
  o="$WORK/rna_$p.bam"
  if [[ ! -s "$o" ]]; then
    echo "RNA total $p: $b"
    # Se reescribe la cabecera para que solo contenga el ADNmt: samtools mpileup empareja las
    # secuencias por su POSICION en la cabecera, y los BAM de ARN pequeno solo tienen el ADNmt.
    # Los mates que cayeron fuera del ADNmt (NUMTs nucleares) se marcan como sin pareja mapeada.
    if bam_indexado "$b"; then
      samtools view -@ "$THREADS" -h "$b" "$MITO_ID"
    else
      samtools view -@ "$THREADS" -h "$b"          # sin indice: se recorre el BAM entero
    fi | awk -v id="$MITO_ID" -v len="$L" 'BEGIN{OFS="\t"}
          /^@HD/ {print; next}
          /^@SQ/ {next}
          /^@/   {next}
          !hdr   {print "@SQ\tSN:" id "\tLN:" len; hdr = 1}
          $3 == id { if ($7 != "=" && $7 != id) {$7 = "*"; $8 = 0; $9 = 0}; print }' \
       | samtools sort -@ "$THREADS" -o "$o" -
    samtools index "$o"
  fi
  n=$(samtools view -c -F 0x904 "$o")
  printf '%s\trnaseq\t%s\tNA\t%s\n' "$p" "$b" "$n" >> "$OUT/conteos.tsv"
  printf '%s\trnaseq\n' "$p" >> "$WORK/nombres.txt"
  echo "  $n lecturas en el ADNmt"
done

# --------------------------------------------------------------- 3. ARN pequeno contra el ADNmt
activate_env "$ENV_SRNA"
if [[ ! -s "$BT1_MITO.1.ebwt" ]]; then
  mkdir -p "$(dirname "$BT1_MITO")"
  bowtie-build -q "$MITO_FA" "$BT1_MITO"
fi
for l in $SRNA_LIBS; do
  f=$(find_srna "$l"); [[ -n "$f" ]] || { echo "ERROR: no encuentro $l en $SRNA_ROOT"; exit 1; }
  o="$WORK/srna_$l.bam"
  if [[ ! -s "$o" ]]; then
    echo "ARN pequeno $l: $f ($(formato_fa "$f"))"
    lecturas_18_35 "$f" > "$WORK/$l.18_35.fa"
    # -v 1: hasta 1 mismatch (las variantes deben poder mapear); -k 1 --best: una posicion por lectura
    bowtie -p "$THREADS" -f -v 1 -k 1 --best --no-unal -S "$BT1_MITO" "$WORK/$l.18_35.fa" 2> "$WORK/$l.bowtie.log" \
      | samtools sort -@ "$THREADS" -o "$o" -
    samtools index "$o"
    n_in=$(grep -c '^>' "$WORK/$l.18_35.fa"); rm -f "$WORK/$l.18_35.fa"
  else
    n_in=NA
  fi
  n=$(samtools view -c -F 0x904 "$o")
  printf '%s\tsrna\t%s\t%s\t%s\n' "$l" "$f" "$n_in" "$n" >> "$OUT/conteos.tsv"
  printf '%s\tsrna\n' "$l" >> "$WORK/nombres.txt"
  echo "  $n de $n_in lecturas de 18-35 nt en el ADNmt"
done

# --------------------------------------------------------------- 4. frecuencias y emparejamiento
bams=()
while read -r n t; do
  [[ "$t" == rnaseq ]] && bams+=("$WORK/rna_$n.bam") || bams+=("$WORK/srna_$n.bam")
done < "$WORK/nombres.txt"
echo
echo "mpileup de ${#bams[@]} BAM sobre $MITO_ID"
# -A: incluir pares "anomalos" (un mate en un NUMT nuclear); -B: sin BAQ; -Q 20: calidad de base
samtools mpileup -A -B -Q 20 -q 0 -d 10000000 -r "$MITO_ID" -f "$MITO_FA" "${bams[@]}" 2> "$WORK/mpileup.log" \
  | python3 "$PY" "$WORK/nombres.txt" "$OUT" --parejas "$PAREJAS" | tee "$OUT/RESUMEN.txt"
if grep -q '\[E::' "$WORK/mpileup.log"; then
  echo "ERROR en mpileup:"; grep '\[E::' "$WORK/mpileup.log" | sort | uniq -c | head; exit 1
fi

cat >> "$OUT/RESUMEN.txt" <<EOF

Lectura: cada libreria de ARN pequeno de control (M2, M11, M21, M60, M72) debe tener en
primer lugar a su propio pool de RNA total. Si "M41" tiene en primer lugar a M48_S5 con una
distancia del mismo orden que las parejas de control, el intercambio queda confirmado.
"M48" no tiene pareja en esta comparacion (M41_S4 fallo); debe quedar lejos de todos.
EOF
echo
echo "Done: $OUT ($((SECONDS/60)) min)"
