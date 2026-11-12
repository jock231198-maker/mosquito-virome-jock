#!/bin/bash
# ---------------------------------------------------------------------------
# srna_map_worker.sh — FASE 2 / NIVEL B, paso 1: particion de las lecturas 18-35
#
# Una libreria por tarea del array. Para cada libreria:
#
#   18-35 ──► rRNA (-v 1) ──► fuera, se cuenta
#            └► resto = DENOMINADOR de RPM ("18-35 no rRNA")
#                 ├─► viral_targets  (-v 1 -a --best --strata)  -> SAM
#                 └─► AaegL5 (-v 0 -k 50 --best --strata), SOLO las lecturas
#                     que dieron hit viral -> tabla lectura/n_loci
#
# POR QUE EL GENOMA VA DESPUES DEL VIRAL Y NO ANTES
#   El plan (PLAN_smallRNA.md §4a) exige que viral y genoma NO vayan en serie:
#   si se quitaran primero las lecturas que mapean a AaegL5 desaparecerian las
#   de los EVE, que son las que se quieren ver. Eso se respeta: NADA se filtra
#   por genoma antes del mapeo viral. Lo unico que cambia aqui es QUE lecturas
#   se preguntan al genoma: solo las que ya tienen hit viral, porque
#   `tambien_en_genoma` solo se usa en esas. Resultado identico, 20-50x mas
#   barato que mapear las 35-45 M de lecturas de cada libreria contra AaegL5.
#
# POR QUE -k 50 Y NO "-a --best --strata -M 50"
#   En bowtie1 -a, -k y -M son alternativas excluyentes: pedir -a y -M a la vez
#   es contradictorio (el plan lo escribia asi). Para el genoma no hace falta la
#   lista de loci, solo (1) si la lectura esta en el genoma y (2) cuantos loci
#   tiene, para que el multimapeo de los EVE en regiones repetidas sea visible
#   (aviso 4 del plan). -k 50 --best --strata da las dos cosas; 50 = "muchos".
#
# PESOS: los crudos del proveedor (00_CleanData/*.clean.fa.gz) vienen
#   COLAPSADOS (>seq_x1234 o ">seq 1234"). No se expanden: cada lectura lleva su
#   multiplicador en el nombre (>r<i>_x<w>) y se arrastra hasta las metricas.
#   Expandir 40 M de lecturas a 750 M de entradas FASTA era lo que hacia B0b
#   (ahi daba igual, eran 7 librerias y un indice minusculo).
#
# Uso:
#   source "$S/srna_env.sh"
#   "$S/srna_map_worker.sh" check                      # en el login
#   N=$(wc -l < "$SCRATCH/nivelB_libs.txt")
#   bsub -J "srnamap[1-$N]%8" \
#        -o "$LOGS_DIR/srnamap.%J.%I.log" -e "$LOGS_DIR/srnamap.%J.%I.err" \
#        -q normal -n 4 -M 8000 -R "select[mem>8000] rusage[mem=8000] span[hosts=1]" \
#        "$S/srna_map_worker.sh run $SCRATCH/nivelB_libs.txt"
#   (recursos = estimacion; medir con bhist -l y anotar en bsubs_reproducibilidad.md)
#
# Variables (exportar antes):
#   BT1_VIRAL        indice viral. Para el catalogo de 500 bp:
#                    export BT1_VIRAL=$SRNA_REFS/bt1_viral_min500/viral_targets
#   SRNA_READS_SRC   proveedor (00_CleanData, por defecto) | propio (trimmed/ de la fase 1)
#   SRNA_LIBS_FILE   lista de librerias, una por linea (por defecto se deduce)
#   KEEP_FA          1 = no borrar los FASTA intermedios de $SCRATCH
#
# SALIDA  $RESULTS_DIR/nivelB/
#   particion/<lib>.particion.tsv     conteos de la particion (lecturas y pesos)
#   metricas/<lib>.votu_metrics.tsv   una fila por vOTU con lecturas en esa libreria
#   sam/<lib>.viral.sam.gz            SAM viral (se conserva: rehacerlo cuesta)
# ---------------------------------------------------------------------------
set -euo pipefail
SECONDS=0
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export S="${S:-$SCRIPT_DIR}"
source "$SCRIPT_DIR/srna_env.sh" "${SRNA_LOTE:-MERIDA_SRNA}" >/dev/null

cmd="${1:?uso: srna_map_worker.sh check|run [lista_librerias]}"
LISTA="${2:-${SRNA_LIBS_FILE:-$SCRATCH/nivelB_libs.txt}}"
THREADS="${LSB_DJOB_NUMPROC:-2}"
SRNA_READS_SRC="${SRNA_READS_SRC:-proveedor}"
KEEP_FA="${KEEP_FA:-0}"

OUT="$RESULTS_DIR/nivelB"
WORK="$SCRATCH/nivelB"

# --------------------------------------------------------------- entradas de lecturas
# proveedor: <SRNA_DATA>/*/00_CleanData/<lib>.clean.fa.gz   (FASTA, colapsado)
# propio   : <SCRATCH>/trimmed/<lib>.fastq.gz               (FASTQ, fase 1)
find_reads() {
  local lib="$1"
  if [[ "$SRNA_READS_SRC" == propio ]]; then
    local f="$SCRATCH/trimmed/$lib.fastq.gz"; [[ -s "$f" ]] && echo "$f" || true
  else
    ls "$SRNA_DATA"/*/00_CleanData/"$lib".clean.fa.gz 2>/dev/null | head -1 || true
  fi
}

formato_fa() {   # colapsado_x | colapsado_espacio | expandido
  local h; h=$( { zcat "$1" 2>/dev/null || true; } | head -1)
  if   [[ "$h" =~ _x[0-9]+$ ]];          then echo "colapsado_x"
  elif [[ "$h" =~ [[:space:]][0-9]+$ ]]; then echo "colapsado_espacio"
  else echo "expandido"; fi
}

# FASTA 18-35 con el peso en el nombre: >r<i>_x<w>.  Imprime a stdout;
# los totales (entradas y peso) van al fichero que se pasa como $2.
extraer_18_35() {
  local f="$1" tot="$2"
  if [[ "$f" == *.fastq.gz || "$f" == *.fq.gz ]]; then
    zcat "$f" | awk -v OUT="$tot" '
      NR % 4 == 2 { L = length($0)
                    if (L >= 18 && L <= 35) { k++; w++; print ">r" k "_x1"; print $0 } }
      END { printf "entradas\t%d\npeso\t%d\n", k, w > OUT }'
  else
    local modo; modo=$(formato_fa "$f")
    zcat "$f" | awk -v modo="$modo" -v OUT="$tot" '
      /^>/ { n = 1
             if (modo == "colapsado_x")            { if (match($0, /_x[0-9]+$/)) n = substr($0, RSTART + 2) + 0 }
             else if (modo == "colapsado_espacio") { n = $NF + 0 }
             if (n < 1) n = 1
             next }
      { L = length($0)
        if (L >= 18 && L <= 35) { k++; w += n; print ">r" k "_x" n; print $0 } }
      END { printf "entradas\t%d\npeso\t%d\n", k, w > OUT }'
  fi
}

peso_fasta() {   # suma de los _x<w> de un FASTA
  awk '/^>/ { n = 1; if (match($0, /_x[0-9]+$/)) n = substr($0, RSTART + 2) + 0; w += n }
       END { print w + 0 }' "$1"
}

# --------------------------------------------------------------- check
if [[ "$cmd" == check ]]; then
  echo "== Nivel B paso 1: comprobacion previa"
  echo
  activate_env "$ENV_SRNA" >/dev/null 2>&1 || echo "  !! no se activa el entorno '$ENV_SRNA'"
  for t in bowtie samtools python3; do
    printf '  %-10s %s\n' "$t" "$(command -v "$t" >/dev/null && "$t" --version 2>&1 | head -1 || echo 'NO ESTA')"
  done
  [[ -s "$SCRIPT_DIR/srna_votu_metrics.py" ]] && echo "  srna_votu_metrics.py ok" \
    || echo "  !! falta $SCRIPT_DIR/srna_votu_metrics.py"
  echo
  echo "-- indices bowtie1"
  for p in "$BT1_RRNA" "$BT1_VIRAL" "$BT1_GENOME"; do
    if [[ -s "$p.1.ebwt" || -s "$p.1.ebwtl" ]]; then
      n=$(bowtie-inspect -n "$p" 2>/dev/null | wc -l || echo '?')
      printf '  [ok ] %-55s %s secuencias\n' "$p" "$n"
      [[ -s "$(dirname "$p")/MANIFEST.txt" ]] && sed 's/^/        /' "$(dirname "$p")/MANIFEST.txt" | head -3
    else
      printf '  [ -- ] %s\n' "$p"
    fi
  done
  echo
  echo "   OJO: el nivel B esta anclado al catalogo de 5,899 vOTUs (corte >=500 bp)."
  echo "   El indice viral se construye con:"
  echo "     export VOTUS_FA=\$REFS_DIR/votu_catalog/MERIDA_cdhit95_85_G1g1_min500/votus.fa"
  echo "     export VIRAL_EXTRA=\$SRNA_REFS/fasta/AEFE1.fa"
  echo "     export BT1_VIRAL=\$SRNA_REFS/bt1_viral_min500/viral_targets"
  echo "     \"\$S/srna_build_refs.sh\" viral       # en bt1_viral_min500, NO pisa el de 1 kb"
  echo
  echo "-- lecturas ($SRNA_READS_SRC)"
  if [[ -s "$LISTA" ]]; then
    n_ok=0; n_no=0
    while read -r lib; do
      [[ -n "$lib" ]] || continue
      f=$(find_reads "$lib")
      if [[ -n "$f" ]]; then
        n_ok=$((n_ok + 1))
        printf '  %-6s %s (%s)\n' "$lib" "${f#"$SRNA_DATA"/}" "$(formato_fa "$f")"
      else
        n_no=$((n_no + 1)); printf '  %-6s !! NO ENCONTRADO\n' "$lib"
      fi
    done < "$LISTA"
    echo "  -> $n_ok encontradas, $n_no sin fichero"
  else
    echo "  !! no existe la lista $LISTA"
    echo "     crearla con los prefijos de libreria, uno por linea; p.ej. desde el samplesheet:"
    echo "       awk -F'\\t' 'NR>1 {print \$1}' \"\$SRNA_SHEET\" > $SCRATCH/nivelB_libs.txt"
  fi
  echo
  echo "-- AEFE1 (control positivo de perfil EVE)"
  if [[ -s "$BT1_VIRAL.1.ebwt" || -s "$BT1_VIRAL.1.ebwtl" ]]; then
    bowtie-inspect -n "$BT1_VIRAL" 2>/dev/null | grep -i -q aefe1 \
      && echo "  ok: AEFE1 esta en el indice viral" \
      || echo "  !! AEFE1 NO esta en el indice: sin control positivo no se puede calibrar la clase EVE"
  else
    echo "  (sin indice viral todavia)"
  fi
  exit 0
fi
[[ "$cmd" == run ]] || { echo "subcomando desconocido: $cmd"; exit 2; }

# --------------------------------------------------------------- run
[[ -s "$LISTA" ]] || { echo "ERROR: no existe la lista $LISTA"; exit 1; }
idx="${LSB_JOBINDEX:-1}"
lib=$(sed -n "${idx}p" "$LISTA")
[[ -n "$lib" ]] || { echo "ERROR: la linea $idx de $LISTA esta vacia"; exit 1; }

for p in "$BT1_RRNA" "$BT1_VIRAL" "$BT1_GENOME"; do
  [[ -s "$p.1.ebwt" || -s "$p.1.ebwtl" ]] || { echo "ERROR: falta el indice $p (corre 'check')"; exit 1; }
done

mkdir -p "$OUT"/{particion,metricas,sam} "$WORK"
echo "nivelB paso 1 — $lib — $(date '+%F %T') en $(hostname), $THREADS hilos"
activate_env "$ENV_SRNA"

f=$(find_reads "$lib")
[[ -n "$f" ]] || { echo "ERROR: no encuentro las lecturas de $lib ($SRNA_READS_SRC)"; exit 1; }
echo "lecturas: $f"

W="$WORK/$lib"
mkdir -p "$W"

# --- 1. 18-35 nt ------------------------------------------------------------
extraer_18_35 "$f" "$W/tot_18_35.tsv" > "$W/in.fa"
n_entradas=$(awk -F'\t' '$1 == "entradas" {print $2}' "$W/tot_18_35.tsv")
n_18_35=$(awk -F'\t' '$1 == "peso" {print $2}' "$W/tot_18_35.tsv")
echo "18-35 nt: $n_18_35 lecturas en $n_entradas secuencias unicas"
(( n_18_35 > 0 )) || { echo "ERROR: ninguna lectura de 18-35 nt en $lib"; exit 1; }

# --- 2. rRNA fuera; el resto es el denominador de RPM -----------------------
# sin --no-unal: en bowtie1 es una opcion de SAM y aqui la salida no es SAM
# (solo interesa --un, el FASTA de lo que NO es rRNA)
bowtie -p "$THREADS" -f -v 1 -k 1 --un "$W/norrna.fa" \
       "$BT1_RRNA" "$W/in.fa" > /dev/null 2> "$W/rrna.log"
n_norrna=$(peso_fasta "$W/norrna.fa")
n_rrna=$((n_18_35 - n_norrna))
printf '18-35 no rRNA: %d (rRNA fuera: %d, %.1f%%)\n' \
  "$n_norrna" "$n_rrna" "$(awk -v a="$n_rrna" -v b="$n_18_35" 'BEGIN{print 100*a/b}')"
(( n_norrna > 0 )) || { echo "ERROR: todo era rRNA en $lib"; exit 1; }

# --- 3. viral: TODAS las lecturas no rRNA, sin filtrar por genoma -----------
bowtie -p "$THREADS" -f -v 1 -a --best --strata -S --no-unal \
       "$BT1_VIRAL" "$W/norrna.fa" 2> "$W/viral.log" | gzip -1 > "$OUT/sam/$lib.viral.sam.gz"
zcat "$OUT/sam/$lib.viral.sam.gz" | awk '!/^@/ {print $1}' | sort -u > "$W/viral_hits.txt"
n_hit=$(awk '{ n = 1; if (match($0, /_x[0-9]+$/)) n = substr($0, RSTART + 2) + 0; w += n }
             END { print w + 0 }' "$W/viral_hits.txt")
echo "con hit viral: $n_hit lecturas en $(wc -l < "$W/viral_hits.txt") secuencias unicas"

# --- 4. genoma: solo las lecturas con hit viral -----------------------------
# tabla lectura -> n_loci (0 = no esta en el genoma). -k 50: 50 significa "muchos".
if [[ -s "$W/viral_hits.txt" ]]; then
  awk 'NR == FNR { q[$1] = 1; next }
       /^>/ { keep = (substr($1, 2) in q); if (keep) print; next }
       keep' "$W/viral_hits.txt" "$W/norrna.fa" > "$W/viral_hits.fa"
  bowtie -p "$THREADS" -f -v 0 -k 50 --best --strata -S --no-unal \
         "$BT1_GENOME" "$W/viral_hits.fa" 2> "$W/genoma.log" \
    | awk -F'\t' '!/^@/ {n[$1]++} END {for (r in n) print r "\t" n[r]}' > "$W/genoma.tsv"
else
  : > "$W/genoma.tsv"
fi
n_gen=$(awk '{ n = 1; if (match($1, /_x[0-9]+$/)) n = substr($1, RSTART + 2) + 0; w += n } END {print w + 0}' "$W/genoma.tsv")
printf 'de esas, tambien en AaegL5: %d lecturas (%.1f%%)\n' \
  "$n_gen" "$(awk -v a="$n_gen" -v b="$n_hit" 'BEGIN{print (b ? 100*a/b : 0)}')"

# --- 5. particion (se conserva: es el denominador y el QC de la libreria) ---
{ printf 'libreria\tmetrica\tvalor\n'
  printf '%s\tlecturas_18_35\t%d\n'        "$lib" "$n_18_35"
  printf '%s\tsecuencias_unicas\t%d\n'     "$lib" "$n_entradas"
  printf '%s\trrna\t%d\n'                  "$lib" "$n_rrna"
  printf '%s\tdenominador_rpm\t%d\n'       "$lib" "$n_norrna"
  printf '%s\tcon_hit_viral\t%d\n'         "$lib" "$n_hit"
  printf '%s\tviral_en_genoma\t%d\n'       "$lib" "$n_gen"
} > "$OUT/particion/$lib.particion.tsv"

# --- 6. metricas por vOTU ---------------------------------------------------
# El join SAM x genoma se hace EN DISCO, no en un diccionario: las lecturas con
# hit viral son millones de nombres distintos y el catalogo de 500 bp es 87 %
# huesped, asi que el lado viral no es pequeno. De paso el flujo queda ordenado
# por nombre de lectura, con lo que "lectura unica" (una sola vOTU) es exacto.
export LC_ALL=C          # sort y join tienen que coincidir en la colacion
longitudes="$W/longitudes.tsv"
zcat "$OUT/sam/$lib.viral.sam.gz" \
  | awk -F'\t' '/^@SQ/ { sn = ""; ln = ""
                         for (i = 2; i <= NF; i++) { if ($i ~ /^SN:/) sn = substr($i, 4)
                                                     if ($i ~ /^LN:/) ln = substr($i, 4) }
                         if (sn != "") print sn "\t" ln; next }
                !/^@/ { exit }' > "$longitudes"

zcat "$OUT/sam/$lib.viral.sam.gz" \
  | awk -F'\t' '!/^@/ { h = (int($2 / 16) % 2) ? "-" : "+"
                        print $1 "\t" $3 "\t" $4 "\t" h "\t" $10 }' \
  | sort -k1,1 -S 2G -T "$W" > "$W/aln.tsv"
sort -k1,1 -S 512M -T "$W" "$W/genoma.tsv" > "$W/genoma.sorted.tsv"

join -t $'\t' -a 1 -e 0 -o '1.1,2.2,1.2,1.3,1.4,1.5' \
     "$W/aln.tsv" "$W/genoma.sorted.tsv" \
  | python3 "$SCRIPT_DIR/srna_votu_metrics.py" \
      --muestra "$lib" --denominador "$n_norrna" --longitudes "$longitudes" \
      --perfil "$OUT/metricas/$lib.perfil_tamanos.tsv" \
      --si "$SRNA_SI_MIN" "$SRNA_SI_MAX" --pi "$SRNA_PI_MIN" "$SRNA_PI_MAX" \
  > "$OUT/metricas/$lib.votu_metrics.tsv"
echo "vOTUs con lecturas: $(( $(wc -l < "$OUT/metricas/$lib.votu_metrics.tsv") - 1 ))"

[[ "$KEEP_FA" == 1 ]] || rm -f "$W"/in.fa "$W"/norrna.fa "$W"/viral_hits.fa \
                               "$W"/aln.tsv "$W"/genoma.sorted.tsv
echo "Done: $lib ($((SECONDS/60)) min)"
